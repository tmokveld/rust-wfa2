use crate::wfa2;
use std::ffi::CString;
use std::io;
use std::os::raw::{c_char, c_void};
use std::panic;
use std::path::Path;

use super::attributes::WFAttributes;
use super::cigar::CigarView;
use super::config::{
    assert_memory_model_compatibility, validate_heuristics_for_distance_metric,
    validate_max_alignment_steps, validate_max_memory, validate_max_num_threads,
    validate_memory_model_compatibility, validate_min_offsets_per_thread, validate_penalties,
    AdaptiveHeuristic, AlignmentResult, AlignmentScope, AlignmentStatus, BandHeuristic,
    DistanceMetric, DropHeuristic, Heuristics, MemoryModel, Penalties, ResourceLimits, WfaError,
};
use super::lambda::{lambda_match_trampoline, LambdaMatcherContext};
use super::packed2bits::validate_packed2bits_sequence;
use super::span::{alignment_span_from_ops, extension_alignment_span_from_ops};

// WFA2 defines DPMATRIX_DIAGONAL_NULL as INT_MAX. Bindgen does not emit
// that macro consistently across libclang/platform combinations.
const DPMATRIX_DIAGONAL_NULL: i32 = i32::MAX;

fn penalties_have_negative_match(penalties: Penalties) -> bool {
    match penalties {
        Penalties::Linear { match_, .. }
        | Penalties::Affine { match_, .. }
        | Penalties::Affine2p { match_, .. } => match_ < 0,
        Penalties::Indel | Penalties::Edit => false,
    }
}

fn path_to_cstring(path: &Path) -> io::Result<CString> {
    use std::os::unix::ffi::OsStrExt;

    CString::new(path.as_os_str().as_bytes()).map_err(|_| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            "plot path contains an interior NUL byte",
        )
    })
}

#[cfg(not(unix))]
fn path_to_cstring(path: &Path) -> io::Result<CString> {
    CString::new(path.as_os_str().to_string_lossy().as_bytes()).map_err(|_| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            "plot path contains an interior NUL byte",
        )
    })
}

pub(crate) struct WfaRawHandle {
    attributes: WFAttributes,
    inner: *mut wfa2::wavefront_aligner_t,
    // Lengths of the last aligned pattern/text. This is the only reliable source for
    // BiWFA (MemoryUltraLow), where the C aligner rewrites its `sequences` bounds during
    // recursion and never restores the originals.
    last_sequence_lengths: Option<(usize, usize)>,
}

impl WfaRawHandle {
    pub(crate) fn new(mut attributes: WFAttributes) -> Result<Self, WfaError> {
        validate_penalties(attributes.penalties())?;
        validate_memory_model_compatibility(
            attributes.selected_memory_model(),
            attributes.alignment_scope_value(),
            attributes.penalties(),
            &Heuristics::none(),
        )?;

        let inner = unsafe { wfa2::wavefront_aligner_new(&mut attributes.inner) };
        Ok(Self {
            attributes,
            inner,
            last_sequence_lengths: None,
        })
    }

    pub(crate) fn alignment_scope(&self) -> AlignmentScope {
        AlignmentScope::from(self.attributes.inner.alignment_scope)
    }

    pub(crate) fn distance_metric(&self) -> DistanceMetric {
        self.attributes.distance_metric()
    }

    pub(crate) fn memory_model(&self) -> MemoryModel {
        self.attributes.selected_memory_model()
    }

    pub(crate) fn penalties(&self) -> Penalties {
        self.attributes.penalties()
    }

    pub(crate) fn heuristics(&self) -> Heuristics {
        let h = &self.attributes.inner.heuristic;
        let strategy = h.strategy;
        let known_strategy = wfa2::wf_heuristic_strategy_wf_heuristic_banded_static
            | wfa2::wf_heuristic_strategy_wf_heuristic_banded_adaptive
            | wfa2::wf_heuristic_strategy_wf_heuristic_wfadaptive
            | wfa2::wf_heuristic_strategy_wf_heuristic_wfmash
            | wfa2::wf_heuristic_strategy_wf_heuristic_xdrop
            | wfa2::wf_heuristic_strategy_wf_heuristic_zdrop;
        if strategy & !known_strategy != 0 {
            panic!("Unknown heuristic strategy: {}", strategy);
        }
        if strategy == wfa2::wf_heuristic_strategy_wf_heuristic_none {
            return Heuristics::new(h.steps_between_cutoffs);
        }

        let mut heuristics = Heuristics::new(h.steps_between_cutoffs);
        if strategy & wfa2::wf_heuristic_strategy_wf_heuristic_wfadaptive != 0 {
            heuristics = heuristics.with_adaptive(AdaptiveHeuristic::WfAdaptive {
                min_wavefront_length: h.min_wavefront_length,
                max_distance_threshold: h.max_distance_threshold,
            });
        } else if strategy & wfa2::wf_heuristic_strategy_wf_heuristic_wfmash != 0 {
            heuristics = heuristics.with_adaptive(AdaptiveHeuristic::WfMash {
                min_wavefront_length: h.min_wavefront_length,
                max_distance_threshold: h.max_distance_threshold,
            });
        }

        if strategy & wfa2::wf_heuristic_strategy_wf_heuristic_xdrop != 0 {
            heuristics = heuristics.with_drop(DropHeuristic::XDrop { xdrop: h.xdrop });
        } else if strategy & wfa2::wf_heuristic_strategy_wf_heuristic_zdrop != 0 {
            heuristics = heuristics.with_drop(DropHeuristic::ZDrop { zdrop: h.zdrop });
        }

        if strategy & wfa2::wf_heuristic_strategy_wf_heuristic_banded_static != 0 {
            heuristics = heuristics.with_band(BandHeuristic::Static {
                min_k: h.min_k,
                max_k: h.max_k,
            });
        } else if strategy & wfa2::wf_heuristic_strategy_wf_heuristic_banded_adaptive != 0 {
            heuristics = heuristics.with_band(BandHeuristic::Adaptive {
                min_k: h.min_k,
                max_k: h.max_k,
            });
        }

        heuristics
    }

    pub(crate) fn resource_limits(&self) -> ResourceLimits {
        self.attributes.resource_limits()
    }

    pub(crate) fn set_alignment_end_to_end(&mut self) {
        unsafe {
            wfa2::wavefront_aligner_set_alignment_end_to_end(self.inner);
        }
    }

    pub(crate) fn set_alignment_ends_free(
        &mut self,
        pattern_begin_free: i32,
        pattern_end_free: i32,
        text_begin_free: i32,
        text_end_free: i32,
    ) {
        unsafe {
            wfa2::wavefront_aligner_set_alignment_free_ends(
                self.inner,
                pattern_begin_free,
                pattern_end_free,
                text_begin_free,
                text_end_free,
            );
        }
    }

    pub(crate) fn set_alignment_extension(&mut self) {
        unsafe {
            wfa2::wavefront_aligner_set_alignment_extension(self.inner);
        }
    }

    /// Panic if ends-free alignment with nonzero free ends was requested under
    /// `MemoryUltraLow` with negative match rewards. WFA2's BiWFA path still
    /// exits the process for that combination.
    fn assert_ends_free_supported(&self, free_ends: [i32; 4]) {
        if self.memory_model() == MemoryModel::MemoryUltraLow
            && free_ends.iter().any(|&free_ends| free_ends != 0)
            && penalties_have_negative_match(self.penalties())
        {
            panic!(
                "Ends-free alignment with negative match rewards is not supported with MemoryUltraLow"
            );
        }
    }

    /// Panic if extension alignment was requested under `MemoryUltraLow`, which
    /// WFA2's BiWFA path exits the process for.
    fn assert_extension_supported(&self) {
        if self.memory_model() == MemoryModel::MemoryUltraLow {
            panic!("Extension alignment is not supported with MemoryUltraLow");
        }
    }

    /// Panic if a lambda/custom matcher was requested under `MemorySingletrack`,
    /// which WFA2 rejects by exiting the process.
    fn assert_lambda_supported(&self) {
        if self.memory_model() == MemoryModel::MemorySingletrack {
            panic!("Lambda/custom sequence inputs are not supported with MemorySingletrack");
        }
    }

    /// Return the underlying aligner pointer, panicking if it is null.
    fn checked_inner(&self) -> *mut wfa2::wavefront_aligner_t {
        if self.inner.is_null() {
            panic!("Internal aligner pointer is null");
        }
        self.inner
    }

    pub(crate) fn align_end_to_end(&mut self, pattern: &[u8], text: &[u8]) -> AlignmentResult {
        self.set_alignment_end_to_end();
        self.align(pattern, text)
    }

    pub(crate) fn align_end_to_end_packed2bits(
        &mut self,
        pattern: &[u8],
        pattern_len: usize,
        text: &[u8],
        text_len: usize,
    ) -> AlignmentResult {
        self.set_alignment_end_to_end();
        self.align_packed2bits(pattern, pattern_len, text, text_len)
    }

    pub(crate) fn align_end_to_end_lambda<F>(
        &mut self,
        pattern_len: usize,
        text_len: usize,
        matcher: &F,
    ) -> AlignmentResult
    where
        F: Fn(usize, usize) -> bool + Sync,
    {
        self.set_alignment_end_to_end();
        self.align_lambda(pattern_len, text_len, matcher)
    }

    pub(crate) fn align_ends_free(
        &mut self,
        pattern: &[u8],
        pattern_begin_free: i32,
        pattern_end_free: i32,
        text: &[u8],
        text_begin_free: i32,
        text_end_free: i32,
    ) -> AlignmentResult {
        self.assert_ends_free_supported([
            pattern_begin_free,
            pattern_end_free,
            text_begin_free,
            text_end_free,
        ]);

        self.set_alignment_ends_free(
            pattern_begin_free,
            pattern_end_free,
            text_begin_free,
            text_end_free,
        );
        self.align(pattern, text)
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn align_ends_free_packed2bits(
        &mut self,
        pattern: &[u8],
        pattern_len: usize,
        pattern_begin_free: i32,
        pattern_end_free: i32,
        text: &[u8],
        text_len: usize,
        text_begin_free: i32,
        text_end_free: i32,
    ) -> AlignmentResult {
        self.assert_ends_free_supported([
            pattern_begin_free,
            pattern_end_free,
            text_begin_free,
            text_end_free,
        ]);

        self.set_alignment_ends_free(
            pattern_begin_free,
            pattern_end_free,
            text_begin_free,
            text_end_free,
        );
        self.align_packed2bits(pattern, pattern_len, text, text_len)
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn align_ends_free_lambda<F>(
        &mut self,
        pattern_len: usize,
        pattern_begin_free: i32,
        pattern_end_free: i32,
        text_len: usize,
        text_begin_free: i32,
        text_end_free: i32,
        matcher: &F,
    ) -> AlignmentResult
    where
        F: Fn(usize, usize) -> bool + Sync,
    {
        self.assert_ends_free_supported([
            pattern_begin_free,
            pattern_end_free,
            text_begin_free,
            text_end_free,
        ]);

        self.set_alignment_ends_free(
            pattern_begin_free,
            pattern_end_free,
            text_begin_free,
            text_end_free,
        );
        self.align_lambda(pattern_len, text_len, matcher)
    }

    pub(crate) fn align_extension(&mut self, pattern: &[u8], text: &[u8]) -> AlignmentResult {
        self.assert_extension_supported();
        self.set_alignment_extension();
        self.align(pattern, text)
    }

    pub(crate) fn align_extension_packed2bits(
        &mut self,
        pattern: &[u8],
        pattern_len: usize,
        text: &[u8],
        text_len: usize,
    ) -> AlignmentResult {
        self.assert_extension_supported();
        self.set_alignment_extension();
        self.align_packed2bits(pattern, pattern_len, text, text_len)
    }

    pub(crate) fn align_extension_lambda<F>(
        &mut self,
        pattern_len: usize,
        text_len: usize,
        matcher: &F,
    ) -> AlignmentResult
    where
        F: Fn(usize, usize) -> bool + Sync,
    {
        self.assert_extension_supported();
        self.set_alignment_extension();
        self.align_lambda(pattern_len, text_len, matcher)
    }

    pub(crate) fn reap(&mut self) {
        unsafe {
            wfa2::wavefront_aligner_reap(self.checked_inner());
        }
    }

    pub(crate) fn plotting_enabled(&self) -> bool {
        self.attributes.inner.plot.enabled
    }

    pub(crate) fn align(&mut self, pattern: &[u8], text: &[u8]) -> AlignmentResult {
        self.last_sequence_lengths = Some((pattern.len(), text.len()));
        let raw_status = unsafe {
            wfa2::wavefront_align(
                self.inner,
                pattern.as_ptr() as *const c_char,
                pattern.len() as i32,
                text.as_ptr() as *const c_char,
                text.len() as i32,
            )
        };
        let result = self.alignment_result();
        debug_assert_eq!(result.status, AlignmentStatus::from(raw_status));
        result
    }

    pub(crate) fn align_packed2bits(
        &mut self,
        pattern: &[u8],
        pattern_len: usize,
        text: &[u8],
        text_len: usize,
    ) -> AlignmentResult {
        let pattern_len_i32 = validate_packed2bits_sequence("pattern", pattern, pattern_len);
        let text_len_i32 = validate_packed2bits_sequence("text", text, text_len);

        self.last_sequence_lengths = Some((pattern_len, text_len));
        let raw_status = unsafe {
            wfa2::wavefront_align_packed2bits(
                self.inner,
                pattern.as_ptr(),
                pattern_len_i32,
                text.as_ptr(),
                text_len_i32,
            )
        };
        let result = self.alignment_result();
        debug_assert_eq!(result.status, AlignmentStatus::from(raw_status));
        result
    }

    pub(crate) fn align_lambda<F>(
        &mut self,
        pattern_len: usize,
        text_len: usize,
        matcher: &F,
    ) -> AlignmentResult
    where
        F: Fn(usize, usize) -> bool + Sync,
    {
        self.assert_lambda_supported();
        self.last_sequence_lengths = Some((pattern_len, text_len));

        let context = LambdaMatcherContext::new(matcher);
        let raw_status = unsafe {
            wfa2::wavefront_align_lambda(
                self.inner,
                Some(lambda_match_trampoline::<F>),
                &context as *const LambdaMatcherContext<'_, F> as *mut c_void,
                pattern_len as i32,
                text_len as i32,
            )
        };

        if let Some(payload) = context.take_panic() {
            panic::resume_unwind(payload);
        }

        let result = self.alignment_result();
        debug_assert_eq!(result.status, AlignmentStatus::from(raw_status));
        result
    }

    pub(crate) fn write_plot<P: AsRef<Path>>(&self, path: P) -> io::Result<()> {
        if !self.plotting_enabled() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "WFA2 plotting was not enabled when the aligner was built",
            ));
        }

        if self.last_sequence_lengths.is_none() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "cannot write a WFA2 plot before running an alignment",
            ));
        }

        if self.inner.is_null() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "internal aligner pointer is null",
            ));
        }

        let c_path = path_to_cstring(path.as_ref())?;
        let mode = CString::new("w").expect("static fopen mode contains no NUL byte");
        let stream = unsafe { wfa2::fopen(c_path.as_ptr(), mode.as_ptr()) };
        if stream.is_null() {
            return Err(io::Error::last_os_error());
        }

        unsafe {
            wfa2::wavefront_plot_print(stream, self.inner);
        }

        let close_result = unsafe { wfa2::fclose(stream) };
        if close_result == 0 {
            Ok(())
        } else {
            Err(io::Error::last_os_error())
        }
    }

    pub(crate) fn alignment_result(&self) -> AlignmentResult {
        let inner = self.checked_inner();
        let status = unsafe { &(*inner).align_status };
        AlignmentResult {
            status: AlignmentStatus::from(status.status),
            score: status.score,
            dropped: status.dropped,
            null_steps: status.num_null_steps,
            memory_used: status.memory_used,
        }
    }

    pub(crate) fn alignment_end_position(&self) -> Option<(usize, usize)> {
        let inner = self.checked_inner();
        let end_pos = unsafe { (*inner).alignment_end_pos };
        if end_pos.k == DPMATRIX_DIAGONAL_NULL || end_pos.offset == wfa2::WAVEFRONT_OFFSET_NULL {
            return None;
        }

        let pattern_end = end_pos.offset as i64 - end_pos.k as i64;
        let text_end = end_pos.offset as i64;
        if pattern_end < 0 || text_end < 0 {
            return None;
        }

        Some((pattern_end as usize, text_end as usize))
    }

    pub(crate) fn alignment_span(&self) -> ((usize, usize), (usize, usize)) {
        let cigar = self
            .cigar_view()
            .expect("CIGAR is null, alignment might have failed or scope was Score");
        let status = self.alignment_result().status;
        let (pattern_end, text_end) = cigar
            .end_position()
            .or_else(|| {
                if self.is_global_alignment() && status == AlignmentStatus::StatusAlgCompleted {
                    Some(self.sequence_lengths())
                } else {
                    None
                }
            })
            .or_else(|| match status {
                AlignmentStatus::StatusAlgCompleted | AlignmentStatus::StatusAlgPartial => {
                    self.alignment_end_position()
                }
                AlignmentStatus::StatusMaxStepsReached
                | AlignmentStatus::StatusOOM
                | AlignmentStatus::StatusUnattainable => None,
            })
            .unwrap_or_else(|| panic!("No valid alignment span is available"));

        if self.is_global_alignment() {
            return ((0, pattern_end), (0, text_end));
        }

        // Extension and ends-free/local: derive the span directly from the active CIGAR so it
        // stays consistent with the reported operations. The endpoint computed above still
        // serves as a validity gate (it panics for failed alignments before we get here), but
        // for these modes the CIGAR is the source of truth. Extension is anchored at the origin
        // (so leading/trailing indels extend the span), whereas ends-free/local strip both to
        // expose just the aligned core. Crucially, this keeps a fully-trimmed extension (empty
        // CIGAR) reporting a `(0, 0)` span instead of the stale wavefront end position.
        if self.is_extension_alignment() {
            return extension_alignment_span_from_ops(cigar.active_operation_bytes());
        }

        alignment_span_from_ops(cigar.active_operation_bytes())
    }

    pub(crate) fn score(&self) -> i32 {
        self.cigar_view()
            .expect("CIGAR is null, alignment might have failed")
            .score
    }

    pub(crate) fn clipped_operation_score(
        operation: char,
        op_length: i32,
        penalties: &Penalties,
    ) -> i32 {
        match penalties {
            Penalties::Indel | Penalties::Edit => match operation {
                'M' => 0,
                'X' | 'D' | 'I' => op_length,
                _ => panic!("Invalid operation: {}", operation),
            },
            Penalties::Linear {
                match_,
                mismatch,
                indel,
            } => -match operation {
                'M' => op_length * match_,
                'X' => op_length * mismatch,
                'D' | 'I' => op_length * indel,
                _ => panic!("Invalid operation: {}", operation),
            },
            Penalties::Affine {
                match_,
                mismatch,
                gap_opening,
                gap_extension,
            } => -match operation {
                'M' => op_length * match_,
                'X' => op_length * mismatch,
                'D' | 'I' => gap_opening + gap_extension * op_length,
                _ => panic!("Invalid operation: {}", operation),
            },
            Penalties::Affine2p {
                match_,
                mismatch,
                gap_opening1,
                gap_extension1,
                gap_opening2,
                gap_extension2,
            } => -match operation {
                'M' => op_length * match_,
                'X' => op_length * mismatch,
                'D' | 'I' => {
                    let score1 = gap_opening1 + gap_extension1 * op_length;
                    let score2 = gap_opening2 + gap_extension2 * op_length;
                    std::cmp::min(score1, score2)
                }
                _ => panic!("Invalid operation: {}", operation),
            },
        }
    }

    pub(crate) fn cigar_score_clipped(&self, flank_len: usize) -> i32 {
        let cigar = self.cigar_view().unwrap();
        let operations = cigar.clipped_operations(flank_len);

        let mut operation_iter = operations.iter().map(|&op| op as u8 as char);
        let Some(mut last_op) = operation_iter.next() else {
            return 0;
        };

        let penalties = self.penalties();
        let mut score = 0;
        let mut op_length = 1;

        for cur_op in operation_iter {
            if cur_op != last_op {
                score += Self::clipped_operation_score(last_op, op_length, &penalties);
                op_length = 0;
            }
            last_op = cur_op;
            op_length += 1;
        }

        score += Self::clipped_operation_score(last_op, op_length, &penalties);
        score
    }

    pub(crate) fn cigar_ptr(&self) -> *mut wfa2::cigar_t {
        if self.inner.is_null() {
            std::ptr::null_mut()
        } else {
            unsafe { (*self.inner).cigar }
        }
    }

    pub(crate) fn cigar_view(&self) -> Option<CigarView<'_>> {
        let cigar = unsafe { self.cigar_ptr().as_ref() }?;
        let operations = if cigar.operations.is_null() || cigar.max_operations <= 0 {
            &[]
        } else {
            unsafe { std::slice::from_raw_parts(cigar.operations, cigar.max_operations as usize) }
        };
        Some(CigarView::new(
            cigar.score,
            cigar.begin_offset,
            cigar.end_offset,
            cigar.end_v,
            cigar.end_h,
            operations,
        ))
    }

    pub(crate) fn active_cigar_bytes(&self) -> Option<&[u8]> {
        let cigar = self.cigar_view()?;
        Some(cigar.active_operation_bytes())
    }

    pub(crate) fn sequence_lengths(&self) -> (usize, usize) {
        // Use the lengths captured at `align` time. Reading them back from the C aligner is
        // unreliable for BiWFA (MemoryUltraLow): the top-level `sequences` is never populated
        // and the bialigner's `wf_forward` is rewritten to sub-problem bounds during recursion.
        self.last_sequence_lengths
            .expect("Sequence lengths are unavailable; no alignment has been performed")
    }

    pub(crate) fn is_global_alignment(&self) -> bool {
        let inner = self.checked_inner();
        unsafe { (*inner).alignment_form.span == wfa2::alignment_span_t_alignment_end2end }
    }

    pub(crate) fn is_extension_alignment(&self) -> bool {
        let inner = self.checked_inner();
        unsafe { (*inner).alignment_form.extension }
    }

    pub(crate) fn wfa_packed_cigar(&self, show_mismatches: bool) -> Vec<u32> {
        if self.inner.is_null() {
            panic!("Internal aligner pointer is null");
        }

        unsafe {
            let mut sam_cigar_buffer_ptr: *mut u32 = std::ptr::null_mut();
            let mut sam_cigar_length: i32 = 0;

            wfa2::cigar_get_CIGAR(
                self.cigar_ptr(),
                show_mismatches,
                &mut sam_cigar_buffer_ptr,
                &mut sam_cigar_length,
            );

            if !sam_cigar_buffer_ptr.is_null() && sam_cigar_length > 0 {
                let cigar_buffer_slice =
                    std::slice::from_raw_parts(sam_cigar_buffer_ptr, sam_cigar_length as usize);
                cigar_buffer_slice.to_vec()
            } else {
                Vec::new()
            }
        }
    }

    pub(crate) fn count_matches(&self) -> i32 {
        if self.inner.is_null() {
            panic!("Internal aligner pointer is null");
        }
        let cigar_ptr = self.cigar_ptr();
        if cigar_ptr.is_null() {
            panic!("CIGAR pointer is null, cannot count matches.");
        }
        unsafe { wfa2::cigar_count_matches(cigar_ptr) }
    }

    pub(crate) fn cigar_score(&mut self) -> i32 {
        let metric = self.distance_metric();
        let cigar = self.cigar_ptr();
        unsafe {
            match metric {
                DistanceMetric::Indel | DistanceMetric::Edit => wfa2::cigar_score_edit(cigar),
                DistanceMetric::GapLinear => {
                    wfa2::cigar_score_gap_linear(cigar, &self.attributes.inner.linear_penalties)
                }
                DistanceMetric::GapAffine => {
                    wfa2::cigar_score_gap_affine(cigar, &self.attributes.inner.affine_penalties)
                }
                DistanceMetric::GapAffine2p => {
                    wfa2::cigar_score_gap_affine2p(cigar, &self.attributes.inner.affine2p_penalties)
                }
            }
        }
    }

    pub(crate) fn set_heuristics(&mut self, heuristics: Heuristics) {
        heuristics.validate();
        validate_heuristics_for_distance_metric(&heuristics, self.distance_metric());
        assert_memory_model_compatibility(
            self.memory_model(),
            self.alignment_scope(),
            self.penalties(),
            &heuristics,
        );

        let cached = &mut self.attributes.inner.heuristic;
        cached.strategy = wfa2::wf_heuristic_strategy_wf_heuristic_none;
        cached.steps_between_cutoffs = heuristics.steps_between_cutoffs();

        unsafe {
            wfa2::wavefront_aligner_set_heuristic_none(self.inner);
        }

        if let Some(adaptive) = heuristics.adaptive() {
            match adaptive {
                AdaptiveHeuristic::WfAdaptive {
                    min_wavefront_length,
                    max_distance_threshold,
                } => {
                    cached.strategy |= wfa2::wf_heuristic_strategy_wf_heuristic_wfadaptive;
                    cached.min_wavefront_length = min_wavefront_length;
                    cached.max_distance_threshold = max_distance_threshold;
                    unsafe {
                        wfa2::wavefront_aligner_set_heuristic_wfadaptive(
                            self.inner,
                            min_wavefront_length,
                            max_distance_threshold,
                            heuristics.steps_between_cutoffs(),
                        );
                    }
                }
                AdaptiveHeuristic::WfMash {
                    min_wavefront_length,
                    max_distance_threshold,
                } => {
                    cached.strategy |= wfa2::wf_heuristic_strategy_wf_heuristic_wfmash;
                    cached.min_wavefront_length = min_wavefront_length;
                    cached.max_distance_threshold = max_distance_threshold;
                    unsafe {
                        wfa2::wavefront_aligner_set_heuristic_wfmash(
                            self.inner,
                            min_wavefront_length,
                            max_distance_threshold,
                            heuristics.steps_between_cutoffs(),
                        );
                    }
                }
            }
        }

        if let Some(drop_heuristic) = heuristics.drop_heuristic() {
            match drop_heuristic {
                DropHeuristic::XDrop { xdrop } => {
                    cached.strategy |= wfa2::wf_heuristic_strategy_wf_heuristic_xdrop;
                    cached.xdrop = xdrop;
                    unsafe {
                        wfa2::wavefront_aligner_set_heuristic_xdrop(
                            self.inner,
                            xdrop,
                            heuristics.steps_between_cutoffs(),
                        );
                    }
                }
                DropHeuristic::ZDrop { zdrop } => {
                    cached.strategy |= wfa2::wf_heuristic_strategy_wf_heuristic_zdrop;
                    cached.zdrop = zdrop;
                    unsafe {
                        wfa2::wavefront_aligner_set_heuristic_zdrop(
                            self.inner,
                            zdrop,
                            heuristics.steps_between_cutoffs(),
                        );
                    }
                }
            }
        }

        if let Some(band) = heuristics.band() {
            match band {
                BandHeuristic::Static { min_k, max_k } => {
                    cached.strategy |= wfa2::wf_heuristic_strategy_wf_heuristic_banded_static;
                    cached.min_k = min_k;
                    cached.max_k = max_k;
                    unsafe {
                        wfa2::wavefront_aligner_set_heuristic_banded_static(
                            self.inner, min_k, max_k,
                        );
                    }
                }
                BandHeuristic::Adaptive { min_k, max_k } => {
                    cached.strategy |= wfa2::wf_heuristic_strategy_wf_heuristic_banded_adaptive;
                    cached.min_k = min_k;
                    cached.max_k = max_k;
                    unsafe {
                        wfa2::wavefront_aligner_set_heuristic_banded_adaptive(
                            self.inner,
                            min_k,
                            max_k,
                            heuristics.steps_between_cutoffs(),
                        );
                    }
                }
            }
        }
    }

    pub(crate) fn set_max_alignment_steps(&mut self, max_alignment_steps: i32) {
        validate_max_alignment_steps(max_alignment_steps);
        self.attributes.inner.system.max_alignment_steps = max_alignment_steps;
        unsafe {
            wfa2::wavefront_aligner_set_max_alignment_steps(self.inner, max_alignment_steps);
        }
    }

    pub(crate) fn set_max_memory(&mut self, max_memory_resident: u64, max_memory_abort: u64) {
        validate_max_memory(max_memory_resident, max_memory_abort);
        self.attributes.inner.system.max_memory_resident = max_memory_resident;
        self.attributes.inner.system.max_memory_abort = max_memory_abort;
        unsafe {
            wfa2::wavefront_aligner_set_max_memory(
                self.inner,
                max_memory_resident,
                max_memory_abort,
            );
        }
    }

    pub(crate) fn set_max_num_threads(&mut self, max_num_threads: i32) {
        validate_max_num_threads(max_num_threads);
        self.attributes.inner.system.max_num_threads = max_num_threads;
        unsafe {
            wfa2::wavefront_aligner_set_max_num_threads(self.inner, max_num_threads);
        }
    }

    pub(crate) fn set_min_offsets_per_thread(&mut self, min_offsets_per_thread: i32) {
        validate_min_offsets_per_thread(min_offsets_per_thread);
        self.attributes.inner.system.min_offsets_per_thread = min_offsets_per_thread;
        unsafe {
            wfa2::wavefront_aligner_set_min_offsets_per_thread(self.inner, min_offsets_per_thread);
        }
    }
}

impl Drop for WfaRawHandle {
    fn drop(&mut self) {
        unsafe {
            if !self.inner.is_null() {
                wfa2::wavefront_aligner_delete(self.inner);
            }
        }
    }
}
