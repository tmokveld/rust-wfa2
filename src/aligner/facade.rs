use std::io;
use std::path::Path;

use super::builder::WFAlignerBuilder;
use super::cigar::{swap_indel_ops_in_cigar_bytes, swap_indel_ops_in_packed_cigar, CigarOp};
use super::config::{
    AlignmentResult, AlignmentScope, Heuristics, MemoryModel, Penalties, ResourceLimits, WfaAlign,
    WfaOp,
};
use super::raw::WfaRawHandle;

pub struct WFAligner {
    pub(crate) raw: WfaRawHandle,
}

impl WFAligner {
    /// Create a builder for configuring a WFAligner
    pub fn builder(alignment_scope: AlignmentScope, memory_model: MemoryModel) -> WFAlignerBuilder {
        WFAlignerBuilder::new(alignment_scope, memory_model)
    }

    pub fn get_penalties(&self) -> Penalties {
        self.raw.penalties()
    }

    pub fn get_heuristics(&self) -> Heuristics {
        self.raw.heuristics()
    }

    pub fn get_resource_limits(&self) -> ResourceLimits {
        self.raw.resource_limits()
    }

    /// Panic when an alignment-only operation is requested on a `Score`-scoped
    /// aligner, which never produces a CIGAR. `action` completes the message
    /// "Cannot {action} when AlignmentScope is Score".
    fn ensure_alignment_scope(&self, action: &str) {
        if self.raw.alignment_scope() == AlignmentScope::Score {
            panic!("Cannot {action} when AlignmentScope is Score");
        }
    }
}

impl WFAligner {
    /// Align byte-slice WFA pattern and WFA text inputs end to end.
    ///
    /// Both inputs are aligned globally over their full lengths. The returned
    /// [`AlignmentResult`] reports WFA2's status and wavefront score snapshot;
    /// with [`AlignmentScope::Alignment`], CIGAR, packed CIGAR, SAM-oriented
    /// CIGAR, match count, clipped score, and alignment span can be read from
    /// the same aligner after this call.
    pub fn align_end_to_end(&mut self, pattern: &[u8], text: &[u8]) -> AlignmentResult {
        self.raw.align_end_to_end(pattern, text)
    }

    /// Align packed 2-bit WFA pattern and WFA text inputs end to end.
    ///
    /// `pattern_len` and `text_len` are logical unpacked sequence lengths, not
    /// byte lengths. The packed slices must contain at least `(len + 3) / 4`
    /// bytes in WFA2's A/C/G/T layout. This method does not pack or allocate
    /// internally; use [`pack_dna_2bits`] when you need to pack ASCII DNA.
    pub fn align_end_to_end_packed2bits(
        &mut self,
        pattern: &[u8],
        pattern_len: usize,
        text: &[u8],
        text_len: usize,
    ) -> AlignmentResult {
        self.raw
            .align_end_to_end_packed2bits(pattern, pattern_len, text, text_len)
    }

    /// Align WFA pattern/text coordinate spaces with a Rust matcher closure.
    ///
    /// `pattern_len` and `text_len` define the half-open WFA pattern/text index
    /// spaces. `matcher` receives zero-based `(pattern_pos, text_pos)`
    /// coordinates and returns whether those positions match. The closure is
    /// borrowed only for this call, but it must be `Sync` because WFA2 may call
    /// it from multiple worker threads when native parallelism is enabled.
    ///
    /// `MemorySingletrack` is rejected because WFA2's singletrack path does not
    /// support lambda/custom matchers.
    ///
    /// Panics inside the matcher are caught at the C callback boundary and
    /// resumed after WFA2 returns.
    pub fn align_end_to_end_lambda<F>(
        &mut self,
        pattern_len: usize,
        text_len: usize,
        matcher: F,
    ) -> AlignmentResult
    where
        F: Fn(usize, usize) -> bool + Sync,
    {
        self.raw
            .align_end_to_end_lambda(pattern_len, text_len, &matcher)
    }

    /// Align byte-slice WFA pattern and WFA text inputs with free ends.
    ///
    /// The free-end counts are expressed on the WFA pattern/text axes and are
    /// forwarded to WFA2's ends-free alignment form. With
    /// [`AlignmentScope::Alignment`], the active CIGAR and alignment span can
    /// be read after this call. `MemoryUltraLow` is rejected when any free-end
    /// count is nonzero, matching WFA2's unsupported BiWFA ends-free path.
    pub fn align_ends_free(
        &mut self,
        pattern: &[u8],
        pattern_begin_free: i32,
        pattern_end_free: i32,
        text: &[u8],
        text_begin_free: i32,
        text_end_free: i32,
    ) -> AlignmentResult {
        self.raw.align_ends_free(
            pattern,
            pattern_begin_free,
            pattern_end_free,
            text,
            text_begin_free,
            text_end_free,
        )
    }

    /// Align packed 2-bit WFA pattern and WFA text inputs with free ends.
    ///
    /// Logical lengths and packed-slice layout follow
    /// [`WFAligner::align_end_to_end_packed2bits`]. Free-end counts use the
    /// same WFA pattern/text axes as [`WFAligner::align_ends_free`].
    /// `MemoryUltraLow` is rejected for nonzero free ends, matching byte-slice
    /// ends-free alignment.
    #[allow(clippy::too_many_arguments)]
    pub fn align_ends_free_packed2bits(
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
        self.raw.align_ends_free_packed2bits(
            pattern,
            pattern_len,
            pattern_begin_free,
            pattern_end_free,
            text,
            text_len,
            text_begin_free,
            text_end_free,
        )
    }

    /// Align WFA pattern/text coordinate spaces with free ends.
    ///
    /// Free-end counts use the same WFA pattern/text axes as
    /// [`WFAligner::align_ends_free`] and otherwise follow the same native
    /// handling. `MemoryUltraLow` is rejected for nonzero free ends, matching
    /// byte-slice ends-free alignment. `MemorySingletrack` is rejected because
    /// WFA2's singletrack path does not support lambda/custom matchers.
    #[allow(clippy::too_many_arguments)]
    pub fn align_ends_free_lambda<F>(
        &mut self,
        pattern_len: usize,
        pattern_begin_free: i32,
        pattern_end_free: i32,
        text_len: usize,
        text_begin_free: i32,
        text_end_free: i32,
        matcher: F,
    ) -> AlignmentResult
    where
        F: Fn(usize, usize) -> bool + Sync,
    {
        self.raw.align_ends_free_lambda(
            pattern_len,
            pattern_begin_free,
            pattern_end_free,
            text_len,
            text_begin_free,
            text_end_free,
            &matcher,
        )
    }

    /// Align a right extension of byte-slice WFA pattern and WFA text inputs.
    ///
    /// WFA2 extension mode anchors the alignment at the origin and trims the
    /// active CIGAR to the maximal-scoring prefix. With
    /// [`AlignmentScope::Alignment`], extension alignments can return
    /// [`AlignmentStatus::StatusAlgPartial`], and the span is derived from the
    /// active CIGAR. `MemoryUltraLow` is rejected because WFA2's BiWFA path
    /// exits the process for extension alignments.
    pub fn align_extension(&mut self, pattern: &[u8], text: &[u8]) -> AlignmentResult {
        self.raw.align_extension(pattern, text)
    }

    /// Align a right extension of packed 2-bit WFA pattern and WFA text inputs.
    ///
    /// Logical lengths and packed-slice layout follow
    /// [`WFAligner::align_end_to_end_packed2bits`]. `MemoryUltraLow` is
    /// rejected because WFA2's BiWFA path exits the process for extension
    /// alignments.
    pub fn align_extension_packed2bits(
        &mut self,
        pattern: &[u8],
        pattern_len: usize,
        text: &[u8],
        text_len: usize,
    ) -> AlignmentResult {
        self.raw
            .align_extension_packed2bits(pattern, pattern_len, text, text_len)
    }

    /// Align a right extension over WFA pattern/text coordinate spaces.
    ///
    /// The matcher contract is the same as
    /// [`WFAligner::align_end_to_end_lambda`]. `MemoryUltraLow` is rejected
    /// because WFA2's BiWFA path exits the process for extension alignments.
    /// `MemorySingletrack` is rejected because WFA2's singletrack path does not
    /// support lambda/custom matchers.
    pub fn align_extension_lambda<F>(
        &mut self,
        pattern_len: usize,
        text_len: usize,
        matcher: F,
    ) -> AlignmentResult
    where
        F: Fn(usize, usize) -> bool + Sync,
    {
        self.raw
            .align_extension_lambda(pattern_len, text_len, &matcher)
    }

    /// Reclaim reusable wavefront memory without destroying the aligner.
    ///
    /// This calls WFA2's explicit memory reclamation hook for buffered
    /// wavefront, slab, and backtrace storage. The aligner's configuration is
    /// preserved, and the next alignment reallocates internal buffers as needed.
    ///
    /// Call this after copying out any alignment result, CIGAR, or derived
    /// alignment data you need to keep. It is valid to call before the first
    /// alignment or multiple times in a row.
    pub fn reap(&mut self) {
        self.raw.reap();
    }

    /// Write WFA2's native `.plot` dump for the last alignment to `path`.
    ///
    /// The aligner must have been built with [`WFAlignerBuilder::with_plotting`]
    /// and an alignment must have been run first. This writes WFA2's upstream
    /// text dump; rendering PNGs is left to external tooling such as WFA2's
    /// `wfa.plot.py` script.
    pub fn write_plot<P: AsRef<Path>>(&self, path: P) -> io::Result<()> {
        self.raw.write_plot(path)
    }

    pub fn score(&self) -> i32 {
        self.raw.score()
    }

    pub fn cigar_score_clipped(&self, flank_len: usize) -> i32 {
        self.ensure_alignment_scope("clip");
        self.raw.cigar_score_clipped(flank_len)
    }

    pub fn set_heuristics(&mut self, heuristics: Heuristics) {
        self.raw.set_heuristics(heuristics);
    }

    pub fn set_max_alignment_steps(&mut self, max_alignment_steps: i32) {
        self.raw.set_max_alignment_steps(max_alignment_steps);
    }

    pub fn set_max_memory(&mut self, max_memory_resident: u64, max_memory_abort: u64) {
        self.raw
            .set_max_memory(max_memory_resident, max_memory_abort);
    }

    pub fn set_max_num_threads(&mut self, max_num_threads: i32) {
        self.raw.set_max_num_threads(max_num_threads);
    }

    pub fn set_min_offsets_per_thread(&mut self, min_offsets_per_thread: i32) {
        self.raw.set_min_offsets_per_thread(min_offsets_per_thread);
    }

    pub fn get_alignment(&self) -> WfaAlign {
        self.ensure_alignment_scope("get alignment");

        let cigar = self.raw.cigar_view().unwrap();
        let raw_operations = cigar.active_operation_bytes();

        let mut operations = Vec::with_capacity(raw_operations.len());
        for &op in raw_operations {
            let operation = WfaOp::from_u8(op);
            operations.push(operation);
        }

        let (pattern_len, text_len) = self.raw.sequence_lengths();
        let ((xstart, xend), (ystart, yend)) = self.raw.alignment_span();

        WfaAlign {
            score: cigar.score,
            ystart,
            yend,
            xstart,
            xend,
            ylen: text_len,
            xlen: pattern_len,
            operations,
        }
    }

    pub fn get_alignment_span(&self) -> ((usize, usize), (usize, usize)) {
        self.ensure_alignment_scope("get alignment span");
        self.raw.alignment_span()
    }

    /// Return WFA2's raw CIGAR operation bytes.
    ///
    /// These operations follow WFA2's native orientation: they describe how to
    /// transform the `pattern` argument into the `text` argument.
    pub fn wfa_cigar_bytes(&self) -> Vec<u8> {
        self.ensure_alignment_scope("get WFA CIGAR bytes");

        let cigar_str = self
            .raw
            .active_cigar_bytes()
            .expect("CIGAR is null, alignment might have failed or scope was Score");
        cigar_str.to_vec()
    }

    /// Return raw CIGAR operation bytes in SAM reference-to-query orientation.
    ///
    /// This method assumes the last alignment used `pattern` as query and `text`
    /// as reference. It converts WFA's pattern-to-text orientation by swapping
    /// `I` and `D` operations.
    pub fn sam_cigar_bytes(&self) -> Vec<u8> {
        let mut cigar = self.wfa_cigar_bytes();
        swap_indel_ops_in_cigar_bytes(&mut cigar);
        cigar
    }

    /// Return WFA2's CIGAR encoded in BAM/SAM's packed integer format.
    ///
    /// The integer packing follows the BAM convention (`len << 4 | op_code`),
    /// but the operation stream keeps WFA2's native orientation: it describes
    /// how to transform the `pattern` argument into the `text` argument.
    pub fn wfa_packed_cigar(&self, show_mismatches: bool) -> Vec<u32> {
        self.ensure_alignment_scope("get WFA packed CIGAR");
        self.raw.wfa_packed_cigar(show_mismatches)
    }

    /// Return a SAM-oriented CIGAR encoded in BAM/SAM's packed integer format.
    ///
    /// This method assumes the last alignment used `pattern` as query and `text`
    /// as reference. It converts WFA's pattern-to-text orientation by swapping
    /// packed `I` and `D` op codes.
    pub fn sam_packed_cigar(&self, show_mismatches: bool) -> Vec<u32> {
        self.ensure_alignment_scope("get SAM CIGAR");
        let mut cigar = self.raw.wfa_packed_cigar(show_mismatches);
        swap_indel_ops_in_packed_cigar(&mut cigar);
        cigar
    }

    /// Decode BAM/SAM packed CIGAR integers into `(length, op)` pairs without
    /// changing operation orientation.
    pub fn decode_packed_cigar(packed_cigar: &[u32]) -> Vec<CigarOp> {
        const SAM_CIGAR_LEN_SHIFT: u32 = 4;
        const SAM_CIGAR_OP_MASK: u32 = 0xF;
        packed_cigar
            .iter()
            .map(|&encoded_op| {
                let len = encoded_op >> SAM_CIGAR_LEN_SHIFT; // Length is in the upper 28 bits
                let op_code = encoded_op & SAM_CIGAR_OP_MASK; // Operation code is in the lower 4 bits
                let op_char = match op_code {
                    0 => 'M', // BAM_CMATCH (Alignment match (can be sequence match or mismatch))
                    1 => 'I', // BAM_CINS (Insertion to the reference)
                    2 => 'D', // BAM_CDEL (Deletion from the reference)
                    3 => 'N', // BAM_CREF_SKIP (Skipped region from the reference)
                    4 => 'S', // BAM_CSOFT_CLIP (Soft clipping (clipped sequences present in SEQ))
                    5 => 'H', // BAM_CHARD_CLIP (Hard clipping (clipped sequences NOT present in SEQ))
                    6 => 'P', // BAM_CPAD (Padding (silent deletion from padded reference))
                    7 => '=', // BAM_CEQUAL (Sequence match)
                    8 => 'X', // BAM_CDIFF (Sequence mismatch)
                    _ => '?', // Unknown operation
                };
                (len as usize, op_char)
            })
            .collect()
    }

    /// Return WFA2's packed CIGAR decoded into `(length, op)` pairs.
    ///
    /// The operation orientation is WFA pattern-to-text.
    pub fn wfa_cigar(&self, show_mismatches: bool) -> Vec<CigarOp> {
        Self::decode_packed_cigar(&self.wfa_packed_cigar(show_mismatches))
    }

    /// Return a SAM-oriented CIGAR decoded into `(length, op)` pairs.
    ///
    /// This method assumes the last alignment used `pattern` as query and `text`
    /// as reference.
    pub fn sam_cigar(&self, show_mismatches: bool) -> Vec<CigarOp> {
        Self::decode_packed_cigar(&self.sam_packed_cigar(show_mismatches))
    }

    /// Counts the number of match ('M') operations in the CIGAR string.
    pub fn count_matches(&self) -> i32 {
        self.ensure_alignment_scope("count matches");
        self.raw.count_matches()
    }

    pub fn cigar_score(&mut self) -> i32 {
        self.ensure_alignment_scope("calculate CIGAR score");
        self.raw.cigar_score()
    }

    #[cfg(test)]
    pub(crate) fn cigar_string(&self, flank_len: Option<usize>) -> String {
        let offset = flank_len.unwrap_or(0);
        let mut cstr = String::new();

        let cigar = self.raw.cigar_view().unwrap();
        let operations = cigar.clipped_operation_bytes(offset);

        let Some((&first_op, remaining_operations)) = operations.split_first() else {
            return cstr;
        };
        let mut last_op = first_op;
        let mut last_op_length = 1;

        for &cur_op in remaining_operations {
            if cur_op == last_op {
                last_op_length += 1;
            } else {
                cstr.push_str(&format!("{}", last_op_length));
                cstr.push(last_op as char);
                last_op = cur_op;
                last_op_length = 1;
            }
        }
        cstr.push_str(&format!("{}", last_op_length));
        cstr.push(last_op as char);
        cstr
    }

    #[cfg(test)]
    pub(crate) fn matching(
        &self,
        pattern: &[u8],
        text: &[u8],
        flank_len: Option<usize>,
    ) -> (String, String, String) {
        let offset = flank_len.unwrap_or(0);

        let mut pattern_iter = pattern.iter().peekable();
        let mut text_iter = text.iter().peekable();

        if offset > 0 {
            text_iter.nth(offset - 1);
            pattern_iter.nth(offset - 1);
        }

        let mut pattern_alg = String::new();
        let mut ops_alg = String::new();
        let mut text_alg = String::new();

        let cigar = self.raw.cigar_view().unwrap();
        let operations = cigar.clipped_operation_bytes(offset);

        for &operation in operations {
            match operation as char {
                'M' => {
                    if pattern_iter.peek() != text_iter.peek() {
                        ops_alg.push('X');
                    } else {
                        ops_alg.push('|');
                    }
                    pattern_alg.push(*pattern_iter.next().unwrap() as char);
                    text_alg.push(*text_iter.next().unwrap() as char);
                }
                'X' => {
                    if pattern_iter.peek() != text_iter.peek() {
                        ops_alg.push(' ');
                    } else {
                        ops_alg.push('X');
                    }
                    pattern_alg.push(*pattern_iter.next().unwrap() as char);
                    text_alg.push(*text_iter.next().unwrap() as char);
                }
                'I' => {
                    pattern_alg.push('-');
                    ops_alg.push(' ');
                    text_alg.push(*text_iter.next().unwrap() as char);
                }
                'D' => {
                    pattern_alg.push(*pattern_iter.next().unwrap() as char);
                    ops_alg.push(' ');
                    text_alg.push('-');
                }
                _ => panic!("Unknown cigar operation"),
            }
        }
        (pattern_alg, ops_alg, text_alg)
    }
}
