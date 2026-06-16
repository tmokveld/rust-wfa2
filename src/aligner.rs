use crate::wfa2;
use std::fmt;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MemoryModel {
    MemoryHigh,
    MemoryMed,
    MemoryLow,
    MemoryUltraLow,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WfaOp {
    Match,
    Subst,
    Ins,
    Del,
}

impl WfaOp {
    fn from_u8(op_char: u8) -> Self {
        match op_char {
            b'M' => WfaOp::Match,
            b'X' => WfaOp::Subst,
            b'I' => WfaOp::Ins,
            b'D' => WfaOp::Del,
            _ => panic!("Invalid alignment operation character {}", op_char),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct WfaAlign {
    /// WFA alignment score (cost).
    pub score: i32,
    /// Start position of alignment in the reference (text). 0-based.
    pub ystart: usize,
    /// Start position of alignment in the query (pattern). 0-based.
    pub xstart: usize,
    /// End position of alignment in the reference (text). 0-based, exclusive.
    pub yend: usize,
    /// End position of alignment in the query (pattern). 0-based, exclusive.
    pub xend: usize,
    /// Length of the reference sequence (text) involved in the alignment.
    pub ylen: usize,
    /// Length of the query sequence (pattern) involved in the alignment.
    pub xlen: usize,
    /// Vector of alignment operations.
    pub operations: Vec<WfaOp>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignmentScope {
    Score,
    Alignment,
}

impl From<u32> for AlignmentScope {
    fn from(value: u32) -> Self {
        match value {
            wfa2::alignment_scope_t_compute_score => AlignmentScope::Score,
            wfa2::alignment_scope_t_compute_alignment => AlignmentScope::Alignment,
            _ => panic!("Unknown alignment scope: {}", value),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdaptiveHeuristic {
    WfAdaptive {
        min_wavefront_length: i32,
        max_distance_threshold: i32,
    },
    WfMash {
        min_wavefront_length: i32,
        max_distance_threshold: i32,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DropHeuristic {
    XDrop { xdrop: i32 },
    ZDrop { zdrop: i32 },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BandHeuristic {
    Static { min_k: i32, max_k: i32 },
    Adaptive { min_k: i32, max_k: i32 },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Heuristics {
    steps_between_cutoffs: i32,
    adaptive: Option<AdaptiveHeuristic>,
    drop_heuristic: Option<DropHeuristic>,
    band: Option<BandHeuristic>,
}

impl Default for Heuristics {
    fn default() -> Self {
        Self::none()
    }
}

impl Heuristics {
    pub fn none() -> Self {
        Self {
            steps_between_cutoffs: 1,
            adaptive: None,
            drop_heuristic: None,
            band: None,
        }
    }

    pub fn new(steps_between_cutoffs: i32) -> Self {
        validate_heuristic_steps(steps_between_cutoffs);
        Self {
            steps_between_cutoffs,
            ..Self::none()
        }
    }

    pub fn wfa2_default() -> Self {
        Self::wf_adaptive(1, 10, 50)
    }

    pub fn wf_adaptive(
        steps_between_cutoffs: i32,
        min_wavefront_length: i32,
        max_distance_threshold: i32,
    ) -> Self {
        Self::new(steps_between_cutoffs).with_adaptive(AdaptiveHeuristic::WfAdaptive {
            min_wavefront_length,
            max_distance_threshold,
        })
    }

    pub fn wf_mash(
        steps_between_cutoffs: i32,
        min_wavefront_length: i32,
        max_distance_threshold: i32,
    ) -> Self {
        Self::new(steps_between_cutoffs).with_adaptive(AdaptiveHeuristic::WfMash {
            min_wavefront_length,
            max_distance_threshold,
        })
    }

    pub fn xdrop(steps_between_cutoffs: i32, xdrop: i32) -> Self {
        Self::new(steps_between_cutoffs).with_drop(DropHeuristic::XDrop { xdrop })
    }

    pub fn zdrop(steps_between_cutoffs: i32, zdrop: i32) -> Self {
        Self::new(steps_between_cutoffs).with_drop(DropHeuristic::ZDrop { zdrop })
    }

    pub fn banded_static(min_k: i32, max_k: i32) -> Self {
        Self::none().with_band(BandHeuristic::Static { min_k, max_k })
    }

    pub fn banded_adaptive(steps_between_cutoffs: i32, min_k: i32, max_k: i32) -> Self {
        Self::new(steps_between_cutoffs).with_band(BandHeuristic::Adaptive { min_k, max_k })
    }

    pub fn with_steps_between_cutoffs(mut self, steps_between_cutoffs: i32) -> Self {
        validate_heuristic_steps(steps_between_cutoffs);
        self.steps_between_cutoffs = steps_between_cutoffs;
        self
    }

    pub fn with_adaptive(mut self, adaptive: AdaptiveHeuristic) -> Self {
        validate_adaptive_heuristic(adaptive);
        self.adaptive = Some(adaptive);
        self
    }

    pub fn with_drop(mut self, drop_heuristic: DropHeuristic) -> Self {
        validate_drop_heuristic(drop_heuristic);
        self.drop_heuristic = Some(drop_heuristic);
        self
    }

    pub fn with_band(mut self, band: BandHeuristic) -> Self {
        validate_band_heuristic(band);
        self.band = Some(band);
        self
    }

    pub fn steps_between_cutoffs(&self) -> i32 {
        self.steps_between_cutoffs
    }

    pub fn adaptive(&self) -> Option<AdaptiveHeuristic> {
        self.adaptive
    }

    pub fn drop_heuristic(&self) -> Option<DropHeuristic> {
        self.drop_heuristic
    }

    pub fn band(&self) -> Option<BandHeuristic> {
        self.band
    }

    pub fn is_none(&self) -> bool {
        self.adaptive.is_none() && self.drop_heuristic.is_none() && self.band.is_none()
    }

    fn validate(&self) {
        validate_heuristic_steps(self.steps_between_cutoffs);
        if let Some(adaptive) = self.adaptive {
            validate_adaptive_heuristic(adaptive);
        }
        if let Some(drop_heuristic) = self.drop_heuristic {
            validate_drop_heuristic(drop_heuristic);
        }
        if let Some(band) = self.band {
            validate_band_heuristic(band);
        }
    }
}

/// Resource controls for bounding WFA2 alignment work.
///
/// If no resource setters are used, WFA2's default attribute values are:
///
/// - `max_alignment_steps`: `i32::MAX` (effectively unlimited)
/// - `max_memory_resident`: `u64::MAX` (WFA2's automatic resident-memory sentinel)
/// - `max_memory_abort`: `u64::MAX` (effectively unlimited)
/// - `max_num_threads`: `1`
/// - `min_offsets_per_thread`: `500`
///
/// `max_num_threads` only affects performance when `wfa2-sys` is built with
/// the `openmp` feature. In testing, OpenMP was workload-sensitive and often
/// flat or slower, so `1` is the safest default.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ResourceLimits {
    /// Maximum WFA score steps before aborting with `StatusMaxStepsReached`.
    ///
    /// Default: `i32::MAX` (effectively unlimited).
    pub max_alignment_steps: i32,
    /// Memory threshold at which the aligner reaps buffered wavefront memory.
    ///
    /// Default: `u64::MAX`, used by WFA2 as its automatic resident-memory
    /// sentinel.
    pub max_memory_resident: u64,
    /// Memory threshold at which the aligner aborts with `StatusOOM`.
    ///
    /// Default: `u64::MAX` (effectively unlimited).
    pub max_memory_abort: u64,
    /// Maximum number of worker threads used by WFA2.
    ///
    /// Default: `1`.
    pub max_num_threads: i32,
    /// Minimum wavefront offsets required before WFA2 starts another worker.
    ///
    /// Default: `500`.
    pub min_offsets_per_thread: i32,
}

impl ResourceLimits {
    pub fn new(
        max_alignment_steps: i32,
        max_memory_resident: u64,
        max_memory_abort: u64,
        max_num_threads: i32,
        min_offsets_per_thread: i32,
    ) -> Self {
        validate_max_alignment_steps(max_alignment_steps);
        validate_max_memory(max_memory_resident, max_memory_abort);
        validate_max_num_threads(max_num_threads);
        validate_min_offsets_per_thread(min_offsets_per_thread);

        Self {
            max_alignment_steps,
            max_memory_resident,
            max_memory_abort,
            max_num_threads,
            min_offsets_per_thread,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DistanceMetric {
    Indel,
    Edit,
    GapLinear,
    GapAffine,
    GapAffine2p,
}

#[derive(Debug, PartialEq, Eq)]
pub enum Penalties {
    Indel, // Conceptually: mismatch=inf, indel=1
    Edit,  // Conceptually: mismatch=1, indel=1
    Linear {
        match_: i32,
        mismatch: i32,
        indel: i32,
    },
    Affine {
        match_: i32,
        mismatch: i32,
        gap_opening: i32,
        gap_extension: i32,
    },
    Affine2p {
        match_: i32,
        mismatch: i32,
        gap_opening1: i32,
        gap_extension1: i32,
        gap_opening2: i32,
        gap_extension2: i32,
    },
}

impl From<u32> for DistanceMetric {
    fn from(value: u32) -> Self {
        match value {
            wfa2::distance_metric_t_indel => DistanceMetric::Indel,
            wfa2::distance_metric_t_edit => DistanceMetric::Edit,
            wfa2::distance_metric_t_gap_linear => DistanceMetric::GapLinear,
            wfa2::distance_metric_t_gap_affine => DistanceMetric::GapAffine,
            wfa2::distance_metric_t_gap_affine_2p => DistanceMetric::GapAffine2p,
            _ => panic!("Unknown distance metric: {}", value),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignmentStatus {
    // OK Status (>=0)
    StatusAlgCompleted = wfa2::WF_STATUS_ALG_COMPLETED as isize,
    StatusAlgPartial = wfa2::WF_STATUS_ALG_PARTIAL as isize,
    // FAILED Status (<0)
    StatusMaxStepsReached = wfa2::WF_STATUS_MAX_STEPS_REACHED as isize,
    StatusOOM = wfa2::WF_STATUS_OOM as isize,
    StatusUnattainable = wfa2::WF_STATUS_UNATTAINABLE as isize,
}

/// Allocation-free snapshot of WFA2's alignment status after an alignment run.
///
/// `score` is WFA2's status/current wavefront score. It can differ from
/// `WFAligner::score()`, which reads the final CIGAR score when one is
/// available. `memory_used` is the current memory usage reported by WFA2 in
/// bytes, not a peak-memory measurement.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AlignmentResult {
    pub status: AlignmentStatus,
    pub score: i32,
    /// Whether WFA2 reports that the alignment was heuristically dropped.
    pub dropped: bool,
    /// Number of contiguous null wavefront steps reported by WFA2.
    pub null_steps: i32,
    /// Current memory usage reported by WFA2, in bytes.
    pub memory_used: u64,
}

impl fmt::Display for AlignmentStatus {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            AlignmentStatus::StatusAlgCompleted => write!(f, "StatusAlgCompleted"),
            AlignmentStatus::StatusAlgPartial => write!(f, "StatusAlgPartial"),
            AlignmentStatus::StatusMaxStepsReached => write!(f, "StatusMaxStepsReached"),
            AlignmentStatus::StatusOOM => write!(f, "StatusOOM"),
            AlignmentStatus::StatusUnattainable => write!(f, "StatusUnattainable"),
        }
    }
}

impl From<i32> for AlignmentStatus {
    fn from(value: i32) -> Self {
        match value {
            x if x == wfa2::WF_STATUS_ALG_COMPLETED as i32 => AlignmentStatus::StatusAlgCompleted,
            x if x == wfa2::WF_STATUS_ALG_PARTIAL as i32 => AlignmentStatus::StatusAlgPartial,
            wfa2::WF_STATUS_MAX_STEPS_REACHED => AlignmentStatus::StatusMaxStepsReached,
            wfa2::WF_STATUS_OOM => AlignmentStatus::StatusOOM,
            wfa2::WF_STATUS_UNATTAINABLE => AlignmentStatus::StatusUnattainable,
            _ => panic!("Unknown alignment status: {}", value),
        }
    }
}

fn validate_max_alignment_steps(max_alignment_steps: i32) {
    assert!(
        max_alignment_steps > 0,
        "max_alignment_steps must be positive"
    );
}

fn validate_max_memory(max_memory_resident: u64, max_memory_abort: u64) {
    assert!(
        max_memory_resident <= max_memory_abort,
        "max_memory_resident must be less than or equal to max_memory_abort"
    );
}

fn validate_max_num_threads(max_num_threads: i32) {
    assert!(max_num_threads > 0, "max_num_threads must be positive");
}

fn validate_min_offsets_per_thread(min_offsets_per_thread: i32) {
    assert!(
        min_offsets_per_thread > 0,
        "min_offsets_per_thread must be positive"
    );
}

fn validate_heuristic_steps(steps_between_cutoffs: i32) {
    assert!(
        steps_between_cutoffs > 0,
        "steps_between_cutoffs must be positive"
    );
}

fn validate_adaptive_heuristic(adaptive: AdaptiveHeuristic) {
    let (min_wavefront_length, max_distance_threshold) = match adaptive {
        AdaptiveHeuristic::WfAdaptive {
            min_wavefront_length,
            max_distance_threshold,
        }
        | AdaptiveHeuristic::WfMash {
            min_wavefront_length,
            max_distance_threshold,
        } => (min_wavefront_length, max_distance_threshold),
    };
    assert!(
        min_wavefront_length > 0,
        "min_wavefront_length must be positive"
    );
    assert!(
        max_distance_threshold >= 0,
        "max_distance_threshold must be non-negative"
    );
}

fn validate_drop_heuristic(drop_heuristic: DropHeuristic) {
    match drop_heuristic {
        DropHeuristic::XDrop { xdrop } => {
            assert!(xdrop >= 0, "xdrop must be non-negative");
        }
        DropHeuristic::ZDrop { zdrop } => {
            assert!(zdrop >= 0, "zdrop must be non-negative");
        }
    }
}

fn validate_band_heuristic(band: BandHeuristic) {
    let (min_k, max_k) = match band {
        BandHeuristic::Static { min_k, max_k } | BandHeuristic::Adaptive { min_k, max_k } => {
            (min_k, max_k)
        }
    };
    assert!(min_k <= max_k, "min_k must be less than or equal to max_k");
}

fn validate_heuristics_for_distance_metric(
    heuristics: &Heuristics,
    distance_metric: DistanceMetric,
) {
    if heuristics.drop_heuristic().is_some()
        && matches!(
            distance_metric,
            DistanceMetric::Indel | DistanceMetric::Edit
        )
    {
        panic!("drop heuristics are not compatible with edit or indel distance metrics");
    }
}

#[derive(Debug, Copy, Clone)]
struct WFAttributes {
    inner: wfa2::wavefront_aligner_attr_t,
}

impl WFAttributes {
    fn default() -> Self {
        let mut inner = unsafe { wfa2::wavefront_aligner_attr_default };
        inner.heuristic.strategy = wfa2::wf_heuristic_strategy_wf_heuristic_none;
        Self { inner }
    }

    fn memory_model(mut self, memory_model: MemoryModel) -> Self {
        let memory_mode = match memory_model {
            MemoryModel::MemoryHigh => wfa2::wavefront_memory_t_wavefront_memory_high,
            MemoryModel::MemoryMed => wfa2::wavefront_memory_t_wavefront_memory_med,
            MemoryModel::MemoryLow => wfa2::wavefront_memory_t_wavefront_memory_low,
            MemoryModel::MemoryUltraLow => wfa2::wavefront_memory_t_wavefront_memory_ultralow,
        };
        self.inner.memory_mode = memory_mode;
        self
    }

    fn alignment_scope(mut self, alignment_scope: AlignmentScope) -> Self {
        let alignment_scope = match alignment_scope {
            AlignmentScope::Score => wfa2::alignment_scope_t_compute_score,
            AlignmentScope::Alignment => wfa2::alignment_scope_t_compute_alignment,
        };
        self.inner.alignment_scope = alignment_scope;
        self
    }

    fn resource_limits(&self) -> ResourceLimits {
        let system = &self.inner.system;
        ResourceLimits {
            max_alignment_steps: system.max_alignment_steps,
            max_memory_resident: system.max_memory_resident,
            max_memory_abort: system.max_memory_abort,
            max_num_threads: system.max_num_threads,
            min_offsets_per_thread: system.min_offsets_per_thread,
        }
    }

    fn set_resource_limits(mut self, resource_limits: ResourceLimits) -> Self {
        self = self.max_alignment_steps(resource_limits.max_alignment_steps);
        self = self.max_memory(
            resource_limits.max_memory_resident,
            resource_limits.max_memory_abort,
        );
        self = self.max_num_threads(resource_limits.max_num_threads);
        self.min_offsets_per_thread(resource_limits.min_offsets_per_thread)
    }

    fn max_alignment_steps(mut self, max_alignment_steps: i32) -> Self {
        validate_max_alignment_steps(max_alignment_steps);
        self.inner.system.max_alignment_steps = max_alignment_steps;
        self
    }

    fn max_memory(mut self, max_memory_resident: u64, max_memory_abort: u64) -> Self {
        validate_max_memory(max_memory_resident, max_memory_abort);
        self.inner.system.max_memory_resident = max_memory_resident;
        self.inner.system.max_memory_abort = max_memory_abort;
        self
    }

    fn max_num_threads(mut self, max_num_threads: i32) -> Self {
        validate_max_num_threads(max_num_threads);
        self.inner.system.max_num_threads = max_num_threads;
        self
    }

    fn min_offsets_per_thread(mut self, min_offsets_per_thread: i32) -> Self {
        validate_min_offsets_per_thread(min_offsets_per_thread);
        self.inner.system.min_offsets_per_thread = min_offsets_per_thread;
        self
    }

    fn indel_penalties(mut self) -> Self {
        self.inner.distance_metric = wfa2::distance_metric_t_indel;
        self
    }

    fn edit_penalties(mut self) -> Self {
        self.inner.distance_metric = wfa2::distance_metric_t_edit;
        self
    }

    fn linear_penalties(mut self, match_: i32, mismatch: i32, indel: i32) -> Self {
        self.inner.distance_metric = wfa2::distance_metric_t_gap_linear;
        self.inner.linear_penalties.match_ = match_; // (Penalty representation usually M <= 0)
        self.inner.linear_penalties.mismatch = mismatch; // (Penalty representation usually X > 0)
        self.inner.linear_penalties.indel = indel; // (Penalty representation usually I > 0)
        self
    }

    fn affine_penalties(
        mut self,
        match_: i32,
        mismatch: i32,
        gap_opening: i32,
        gap_extension: i32,
    ) -> Self {
        self.inner.distance_metric = wfa2::distance_metric_t_gap_affine;
        self.inner.affine_penalties.match_ = match_; // (Penalty representation usually M <= 0)
        self.inner.affine_penalties.mismatch = mismatch; // (Penalty representation usually X > 0)
        self.inner.affine_penalties.gap_opening = gap_opening; // (Penalty representation usually O > 0)
        self.inner.affine_penalties.gap_extension = gap_extension; // (Penalty representation usually E > 0)
        self
    }

    fn affine2p_penalties(
        mut self,
        match_: i32,
        mismatch: i32,
        gap_opening1: i32,
        gap_extension1: i32,
        gap_opening2: i32,
        gap_extension2: i32,
    ) -> Self {
        self.inner.distance_metric = wfa2::distance_metric_t_gap_affine_2p;
        self.inner.affine2p_penalties.match_ = match_; // (Penalty representation usually M <= 0)
        self.inner.affine2p_penalties.mismatch = mismatch; // (Penalty representation usually X > 0)
                                                           // Usually concave Q1 + E1 < Q2 + E2 and E1 > E2.
        self.inner.affine2p_penalties.gap_opening1 = gap_opening1; // (Penalty representation usually O1 > 0)
        self.inner.affine2p_penalties.gap_extension1 = gap_extension1; // (Penalty representation usually E1 > 0)
        self.inner.affine2p_penalties.gap_opening2 = gap_opening2; // (Penalty representation usually O2 > 0)
        self.inner.affine2p_penalties.gap_extension2 = gap_extension2; // (Penalty representation usually E2 > 0)
        self
    }
}

pub struct WFAlignerBuilder {
    attributes: WFAttributes,
    penalty_set: bool,
    heuristics: Option<Heuristics>,
}

impl WFAlignerBuilder {
    pub fn new(alignment_scope: AlignmentScope, memory_model: MemoryModel) -> Self {
        let attributes = WFAttributes::default()
            .memory_model(memory_model)
            .alignment_scope(alignment_scope);
        Self {
            attributes,
            penalty_set: false,
            heuristics: None,
        }
    }

    /// Configure for indel penalties (Longest Common Subsequence - LCS)
    pub fn indel(mut self) -> Self {
        self.attributes = self.attributes.indel_penalties();
        self.penalty_set = true;
        self
    }

    /// Configure for edit penalties (Levenshtein)
    pub fn edit(mut self) -> Self {
        self.attributes = self.attributes.edit_penalties();
        self.penalty_set = true;
        self
    }

    /// Configure for gap-linear penalties (Needleman-Wunsch) with match_ = 0
    pub fn linear(self, mismatch: i32, indel: i32) -> Self {
        self.linear_with_match(0, mismatch, indel)
    }

    /// Configure for gap-linear penalties (Needleman-Wunsch) with explicit match score
    pub fn linear_with_match(mut self, match_: i32, mismatch: i32, indel: i32) -> Self {
        self.attributes = self.attributes.linear_penalties(match_, mismatch, indel);
        self.penalty_set = true;
        self
    }

    /// Configure for gap-affine penalties (Smith-Waterman-Gotoh) with match_ = 0
    pub fn affine(self, mismatch: i32, gap_opening: i32, gap_extension: i32) -> Self {
        self.affine_with_match(0, mismatch, gap_opening, gap_extension)
    }

    /// Configure for gap-affine penalties (Smith-Waterman-Gotoh) with explicit match score
    pub fn affine_with_match(
        mut self,
        match_: i32,
        mismatch: i32,
        gap_opening: i32,
        gap_extension: i32,
    ) -> Self {
        self.attributes =
            self.attributes
                .affine_penalties(match_, mismatch, gap_opening, gap_extension);
        self.penalty_set = true;
        self
    }

    /// Configure for gap-affine dual-cost penalties (concave 2-pieces) with match_ = 0
    pub fn affine2p(
        self,
        mismatch: i32,
        gap_opening1: i32,
        gap_extension1: i32,
        gap_opening2: i32,
        gap_extension2: i32,
    ) -> Self {
        self.affine2p_with_match(
            0,
            mismatch,
            gap_opening1,
            gap_extension1,
            gap_opening2,
            gap_extension2,
        )
    }

    /// Configure for gap-affine dual-cost penalties (concave 2-pieces) with explicit match score
    #[allow(clippy::too_many_arguments)]
    pub fn affine2p_with_match(
        mut self,
        match_: i32,
        mismatch: i32,
        gap_opening1: i32,
        gap_extension1: i32,
        gap_opening2: i32,
        gap_extension2: i32,
    ) -> Self {
        self.attributes = self.attributes.affine2p_penalties(
            match_,
            mismatch,
            gap_opening1,
            gap_extension1,
            gap_opening2,
            gap_extension2,
        );
        self.penalty_set = true;
        self
    }

    /// Set heuristic configuration for the aligner.
    pub fn with_heuristics(mut self, heuristics: Heuristics) -> Self {
        heuristics.validate();
        self.heuristics = Some(heuristics);
        self
    }

    /// Set all resource limits at once.
    pub fn with_resource_limits(mut self, resource_limits: ResourceLimits) -> Self {
        self.attributes = self.attributes.set_resource_limits(resource_limits);
        self
    }

    /// Set the maximum WFA score steps before aborting with `StatusMaxStepsReached`.
    pub fn with_max_alignment_steps(mut self, max_alignment_steps: i32) -> Self {
        self.attributes = self.attributes.max_alignment_steps(max_alignment_steps);
        self
    }

    /// Set WFA2 memory thresholds in bytes.
    ///
    /// `max_memory_resident` controls when resident buffered memory is reaped.
    /// `max_memory_abort` controls when alignment aborts with `StatusOOM`.
    pub fn with_max_memory(mut self, max_memory_resident: u64, max_memory_abort: u64) -> Self {
        self.attributes = self
            .attributes
            .max_memory(max_memory_resident, max_memory_abort);
        self
    }

    /// Set the maximum number of worker threads used by WFA2.
    pub fn with_max_num_threads(mut self, max_num_threads: i32) -> Self {
        self.attributes = self.attributes.max_num_threads(max_num_threads);
        self
    }

    /// Set the minimum wavefront offsets required before WFA2 starts another worker.
    pub fn with_min_offsets_per_thread(mut self, min_offsets_per_thread: i32) -> Self {
        self.attributes = self
            .attributes
            .min_offsets_per_thread(min_offsets_per_thread);
        self
    }

    /// Build the WFAligner with the configured settings
    pub fn build(self) -> WFAligner {
        if !self.penalty_set {
            panic!("Must set a penalty model before building the aligner");
        }

        let mut raw = WfaRawHandle::new(self.attributes);

        if let Some(heuristics) = self.heuristics {
            raw.set_heuristics(heuristics);
        }

        WFAligner { raw }
    }
}

struct CigarView<'a> {
    score: i32,
    begin_offset: i32,
    end_offset: i32,
    end_v: i32,
    end_h: i32,
    operations: &'a [std::os::raw::c_char],
}

impl CigarView<'_> {
    fn active_operation_bytes(&self) -> &[u8] {
        let operations = &self.operations[self.begin_offset as usize..self.end_offset as usize];
        // SAFETY: i8/u8 same layout, slice borrowed for 'self
        unsafe { std::slice::from_raw_parts(operations.as_ptr() as *const u8, operations.len()) }
    }

    fn end_position(&self) -> Option<(usize, usize)> {
        if self.end_v < 0 || self.end_h < 0 {
            return None;
        }

        Some((self.end_v as usize, self.end_h as usize))
    }
}

/// Derives the aligned span (start and end on both axes) from the active CIGAR operations.
///
/// Both leading and trailing indels are stripped, so the span covers only the aligned core
/// (from the first to the last `M`/`X` column). This keeps the two axes symmetric and is used
/// for ends-free/local alignments, which are always computed unidirectionally (BiWFA is rejected
/// for those modes), so the full CIGAR is available here.
fn alignment_span_from_ops(raw_operations: &[u8]) -> ((usize, usize), (usize, usize)) {
    let mut pattern_index = 0;
    let mut text_index = 0;

    let mut pattern_start = None;
    let mut text_start = None;
    let mut pattern_end = 0;
    let mut text_end = 0;

    for &op in raw_operations {
        match op {
            b'I' => {
                text_index += 1;
            }
            b'D' => {
                pattern_index += 1;
            }
            b'M' | b'X' => {
                pattern_start.get_or_insert(pattern_index);
                text_start.get_or_insert(text_index);

                pattern_index += 1;
                text_index += 1;

                pattern_end = pattern_index;
                text_end = text_index;
            }
            _ => panic!("unexpected WFA operation: {}", op as char),
        }
    }

    (
        (pattern_start.unwrap_or(0), pattern_end),
        (text_start.unwrap_or(0), text_end),
    )
}

/// Derives the aligned span for an extension alignment from the active CIGAR operations.
///
/// Extension alignments are anchored at the origin, so the span always starts at `(0, 0)` and
/// the end is simply the number of pattern/text characters consumed by the CIGAR. Unlike
/// [`alignment_span_from_ops`], leading and trailing indels are *not* stripped (they advance the
/// end on their axis). Deriving the span from the CIGAR rather than the wavefront end position
/// keeps it consistent with the reported operations even when the maximal-scoring prefix is
/// empty (a fully-trimmed extension yields an empty CIGAR and therefore a `(0, 0)` span).
fn extension_alignment_span_from_ops(raw_operations: &[u8]) -> ((usize, usize), (usize, usize)) {
    let mut pattern_end = 0;
    let mut text_end = 0;

    for &op in raw_operations {
        match op {
            b'I' => {
                text_end += 1;
            }
            b'D' => {
                pattern_end += 1;
            }
            b'M' | b'X' => {
                pattern_end += 1;
                text_end += 1;
            }
            _ => panic!("unexpected WFA operation: {}", op as char),
        }
    }

    ((0, pattern_end), (0, text_end))
}

struct WfaRawHandle {
    attributes: WFAttributes,
    inner: *mut wfa2::wavefront_aligner_t,
    // Lengths of the last aligned pattern/text. This is the only reliable source for
    // BiWFA (MemoryUltraLow), where the C aligner rewrites its `sequences` bounds during
    // recursion and never restores the originals.
    last_sequence_lengths: Option<(usize, usize)>,
}

impl WfaRawHandle {
    fn new(mut attributes: WFAttributes) -> Self {
        let inner = unsafe { wfa2::wavefront_aligner_new(&mut attributes.inner) };
        Self {
            attributes,
            inner,
            last_sequence_lengths: None,
        }
    }

    fn alignment_scope(&self) -> AlignmentScope {
        AlignmentScope::from(self.attributes.inner.alignment_scope)
    }

    fn distance_metric(&self) -> DistanceMetric {
        DistanceMetric::from(self.attributes.inner.distance_metric)
    }

    fn memory_model(&self) -> MemoryModel {
        match self.attributes.inner.memory_mode {
            wfa2::wavefront_memory_t_wavefront_memory_high => MemoryModel::MemoryHigh,
            wfa2::wavefront_memory_t_wavefront_memory_med => MemoryModel::MemoryMed,
            wfa2::wavefront_memory_t_wavefront_memory_low => MemoryModel::MemoryLow,
            wfa2::wavefront_memory_t_wavefront_memory_ultralow => MemoryModel::MemoryUltraLow,
            _ => panic!(
                "Unknown memory model: {}",
                self.attributes.inner.memory_mode
            ),
        }
    }

    fn penalties(&self) -> Penalties {
        match self.distance_metric() {
            DistanceMetric::Indel => Penalties::Indel,
            DistanceMetric::Edit => Penalties::Edit,
            DistanceMetric::GapLinear => Penalties::Linear {
                match_: self.attributes.inner.linear_penalties.match_,
                mismatch: self.attributes.inner.linear_penalties.mismatch,
                indel: self.attributes.inner.linear_penalties.indel,
            },
            DistanceMetric::GapAffine => Penalties::Affine {
                match_: self.attributes.inner.affine_penalties.match_,
                mismatch: self.attributes.inner.affine_penalties.mismatch,
                gap_opening: self.attributes.inner.affine_penalties.gap_opening,
                gap_extension: self.attributes.inner.affine_penalties.gap_extension,
            },
            DistanceMetric::GapAffine2p => Penalties::Affine2p {
                match_: self.attributes.inner.affine2p_penalties.match_,
                mismatch: self.attributes.inner.affine2p_penalties.mismatch,
                gap_opening1: self.attributes.inner.affine2p_penalties.gap_opening1,
                gap_extension1: self.attributes.inner.affine2p_penalties.gap_extension1,
                gap_opening2: self.attributes.inner.affine2p_penalties.gap_opening2,
                gap_extension2: self.attributes.inner.affine2p_penalties.gap_extension2,
            },
        }
    }

    fn heuristics(&self) -> Heuristics {
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

    fn resource_limits(&self) -> ResourceLimits {
        self.attributes.resource_limits()
    }

    fn set_alignment_end_to_end(&mut self) {
        unsafe {
            wfa2::wavefront_aligner_set_alignment_end_to_end(self.inner);
        }
    }

    fn set_alignment_ends_free(
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

    fn set_alignment_extension(&mut self) {
        unsafe {
            wfa2::wavefront_aligner_set_alignment_extension(self.inner);
        }
    }

    fn align_end_to_end(&mut self, pattern: &[u8], text: &[u8]) -> AlignmentResult {
        self.set_alignment_end_to_end();
        self.align(pattern, text)
    }

    fn align_ends_free(
        &mut self,
        pattern: &[u8],
        pattern_begin_free: i32,
        pattern_end_free: i32,
        text: &[u8],
        text_begin_free: i32,
        text_end_free: i32,
    ) -> AlignmentResult {
        if self.memory_model() == MemoryModel::MemoryUltraLow
            && [
                pattern_begin_free,
                pattern_end_free,
                text_begin_free,
                text_end_free,
            ]
            .iter()
            .any(|&free_ends| free_ends != 0)
        {
            panic!("Ends-free alignment is not supported with MemoryUltraLow");
        }

        self.set_alignment_ends_free(
            pattern_begin_free,
            pattern_end_free,
            text_begin_free,
            text_end_free,
        );
        self.align(pattern, text)
    }

    fn align_extension(&mut self, pattern: &[u8], text: &[u8]) -> AlignmentResult {
        if self.memory_model() == MemoryModel::MemoryUltraLow {
            panic!("Extension alignment is not supported with MemoryUltraLow");
        }

        self.set_alignment_extension();
        self.align(pattern, text)
    }

    fn reap(&mut self) {
        if self.inner.is_null() {
            panic!("Internal aligner pointer is null");
        }

        unsafe {
            wfa2::wavefront_aligner_reap(self.inner);
        }
    }

    fn align(&mut self, pattern: &[u8], text: &[u8]) -> AlignmentResult {
        self.last_sequence_lengths = Some((pattern.len(), text.len()));
        let raw_status = unsafe {
            wfa2::wavefront_align(
                self.inner,
                pattern.as_ptr() as *const i8,
                pattern.len() as i32,
                text.as_ptr() as *const i8,
                text.len() as i32,
            )
        };
        let result = self.alignment_result();
        debug_assert_eq!(result.status, AlignmentStatus::from(raw_status));
        result
    }

    fn alignment_result(&self) -> AlignmentResult {
        if self.inner.is_null() {
            panic!("Internal aligner pointer is null");
        }

        let status = unsafe { &(*self.inner).align_status };
        AlignmentResult {
            status: AlignmentStatus::from(status.status),
            score: status.score,
            dropped: status.dropped,
            null_steps: status.num_null_steps,
            memory_used: status.memory_used,
        }
    }

    fn alignment_end_position(&self) -> Option<(usize, usize)> {
        if self.inner.is_null() {
            panic!("Internal aligner pointer is null");
        }

        let end_pos = unsafe { (*self.inner).alignment_end_pos };
        if end_pos.k == wfa2::DPMATRIX_DIAGONAL_NULL as i32
            || end_pos.offset == wfa2::WAVEFRONT_OFFSET_NULL
        {
            return None;
        }

        let pattern_end = end_pos.offset as i64 - end_pos.k as i64;
        let text_end = end_pos.offset as i64;
        if pattern_end < 0 || text_end < 0 {
            return None;
        }

        Some((pattern_end as usize, text_end as usize))
    }

    fn alignment_span(&self) -> ((usize, usize), (usize, usize)) {
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

    fn score(&self) -> i32 {
        self.cigar_view()
            .expect("CIGAR is null, alignment might have failed")
            .score
    }

    fn clipped_operation_score(operation: char, op_length: i32, penalties: &Penalties) -> i32 {
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

    fn cigar_score_clipped(&self, flank_len: usize) -> i32 {
        let cigar = self.cigar_view().unwrap();
        let begin_offset = cigar.begin_offset as isize + flank_len as isize;
        let end_offset =
            std::cmp::max(begin_offset, cigar.end_offset as isize - flank_len as isize);

        if begin_offset >= end_offset {
            return 0;
        }

        let penalties = self.penalties();
        let mut score = 0;
        let mut op_length = 0;
        let mut last_op: Option<char> = Some(cigar.operations[begin_offset as usize] as u8 as char);

        for i in begin_offset..end_offset {
            let cur_op = cigar.operations[i as usize] as u8 as char;
            if Some(cur_op) != last_op {
                if let Some(op) = last_op {
                    score += Self::clipped_operation_score(op, op_length, &penalties);
                }
                op_length = 0;
            }
            last_op = Some(cur_op);
            op_length += 1;
        }

        if let Some(op) = last_op {
            score += Self::clipped_operation_score(op, op_length, &penalties);
        }
        score
    }

    fn cigar_ptr(&self) -> *mut wfa2::cigar_t {
        if self.inner.is_null() {
            std::ptr::null_mut()
        } else {
            unsafe { (*self.inner).cigar }
        }
    }

    fn cigar_view(&self) -> Option<CigarView<'_>> {
        let cigar = unsafe { self.cigar_ptr().as_ref() }?;
        let operations = if cigar.operations.is_null() || cigar.max_operations <= 0 {
            &[]
        } else {
            unsafe { std::slice::from_raw_parts(cigar.operations, cigar.max_operations as usize) }
        };
        Some(CigarView {
            score: cigar.score,
            begin_offset: cigar.begin_offset,
            end_offset: cigar.end_offset,
            end_v: cigar.end_v,
            end_h: cigar.end_h,
            operations,
        })
    }

    fn active_cigar_bytes(&self) -> Option<&[u8]> {
        let cigar = unsafe { self.cigar_ptr().as_ref() }?;
        if cigar.operations.is_null() || cigar.begin_offset > cigar.end_offset {
            return Some(&[]);
        }

        let cigar_length = (cigar.end_offset - cigar.begin_offset) as usize;
        let cigar_operations =
            unsafe { cigar.operations.offset(cigar.begin_offset as isize) as *const u8 };
        Some(unsafe { std::slice::from_raw_parts(cigar_operations, cigar_length) })
    }

    fn sequence_lengths(&self) -> (usize, usize) {
        // Use the lengths captured at `align` time. Reading them back from the C aligner is
        // unreliable for BiWFA (MemoryUltraLow): the top-level `sequences` is never populated
        // and the bialigner's `wf_forward` is rewritten to sub-problem bounds during recursion.
        self.last_sequence_lengths
            .expect("Sequence lengths are unavailable; no alignment has been performed")
    }

    fn is_global_alignment(&self) -> bool {
        if self.inner.is_null() {
            panic!("Internal aligner pointer is null");
        }
        unsafe { (*self.inner).alignment_form.span == wfa2::alignment_span_t_alignment_end2end }
    }

    fn is_extension_alignment(&self) -> bool {
        if self.inner.is_null() {
            panic!("Internal aligner pointer is null");
        }
        unsafe { (*self.inner).alignment_form.extension }
    }

    fn sam_cigar(&self, show_mismatches: bool) -> Vec<u32> {
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

    fn count_matches(&self) -> i32 {
        if self.inner.is_null() {
            panic!("Internal aligner pointer is null");
        }
        let cigar_ptr = self.cigar_ptr();
        if cigar_ptr.is_null() {
            panic!("CIGAR pointer is null, cannot count matches.");
        }
        unsafe { wfa2::cigar_count_matches(cigar_ptr) }
    }

    fn cigar_score(&mut self) -> i32 {
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

    fn set_heuristics(&mut self, heuristics: Heuristics) {
        heuristics.validate();
        validate_heuristics_for_distance_metric(&heuristics, self.distance_metric());

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

    fn set_max_alignment_steps(&mut self, max_alignment_steps: i32) {
        validate_max_alignment_steps(max_alignment_steps);
        self.attributes.inner.system.max_alignment_steps = max_alignment_steps;
        unsafe {
            wfa2::wavefront_aligner_set_max_alignment_steps(self.inner, max_alignment_steps);
        }
    }

    fn set_max_memory(&mut self, max_memory_resident: u64, max_memory_abort: u64) {
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

    fn set_max_num_threads(&mut self, max_num_threads: i32) {
        validate_max_num_threads(max_num_threads);
        self.attributes.inner.system.max_num_threads = max_num_threads;
        unsafe {
            wfa2::wavefront_aligner_set_max_num_threads(self.inner, max_num_threads);
        }
    }

    fn set_min_offsets_per_thread(&mut self, min_offsets_per_thread: i32) {
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

// TODO: Unify different Cigar wrappers
/// Represents a single operation: (length, op).
pub type CigarOp = (usize, char);

pub struct WFAligner {
    raw: WfaRawHandle,
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
}

impl WFAligner {
    pub fn align_end_to_end(&mut self, pattern: &[u8], text: &[u8]) -> AlignmentResult {
        self.raw.align_end_to_end(pattern, text)
    }

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

    /// Align a right extension using WFA2's extension mode.
    ///
    /// With `AlignmentScope::Alignment`, WFA2 trims the active CIGAR to the
    /// maximal-scoring extension and can return `StatusAlgPartial`.
    /// `MemoryUltraLow` is rejected because WFA2's BiWFA path exits the process
    /// for extension alignments.
    pub fn align_extension(&mut self, pattern: &[u8], text: &[u8]) -> AlignmentResult {
        self.raw.align_extension(pattern, text)
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

    pub fn score(&self) -> i32 {
        self.raw.score()
    }

    pub fn cigar_score_clipped(&self, flank_len: usize) -> i32 {
        if self.raw.alignment_scope() == AlignmentScope::Score {
            panic!("Cannot clip when AlignmentScope is Score");
        }
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
        if self.raw.alignment_scope() == AlignmentScope::Score {
            panic!("Cannot get alignment when AlignmentScope is Score");
        }

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
        if self.raw.alignment_scope() == AlignmentScope::Score {
            panic!("Cannot get alignment span when AlignmentScope is Score");
        }

        self.raw.alignment_span()
    }

    pub fn cigar_operations(&self) -> Vec<u8> {
        if self.raw.alignment_scope() == AlignmentScope::Score {
            return Vec::new();
        }

        let cigar_str = self
            .raw
            .active_cigar_bytes()
            .expect("CIGAR is null, alignment might have failed or scope was Score");
        cigar_str.to_vec()
    }

    pub fn get_sam_cigar(&self, show_mismatches: bool) -> Vec<u32> {
        if self.raw.alignment_scope() == AlignmentScope::Score {
            panic!("Cannot get SAM CIGAR when AlignmentScope is Score");
        }
        self.raw.sam_cigar(show_mismatches)
    }

    pub fn decode_sam_cigar(sam_cigar_buffer: &[u32]) -> Vec<CigarOp> {
        const SAM_CIGAR_LEN_SHIFT: u32 = 4;
        const SAM_CIGAR_OP_MASK: u32 = 0xF;
        sam_cigar_buffer
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

    /// Counts the number of match ('M') operations in the CIGAR string.
    pub fn count_matches(&self) -> i32 {
        if self.raw.alignment_scope() == AlignmentScope::Score {
            panic!("Cannot count matches when AlignmentScope is Score");
        }
        self.raw.count_matches()
    }

    pub fn cigar_score(&mut self) -> i32 {
        if self.raw.alignment_scope() == AlignmentScope::Score {
            panic!("Cannot calculate CIGAR score when AlignmentScope is Score");
        }
        self.raw.cigar_score()
    }

    #[allow(dead_code)]
    fn cigar_string(&self, flank_len: Option<usize>) -> String {
        let offset = flank_len.unwrap_or(0);
        let mut cstr = String::new();

        let cigar = self.raw.cigar_view().unwrap();

        let begin_offset = cigar.begin_offset as usize + offset;
        let end_offset = cigar.end_offset as usize - offset;

        if begin_offset >= end_offset {
            return cstr;
        }

        let operations = cigar.operations;
        let mut last_op = operations[begin_offset];
        let mut last_op_length = 1;

        for i in 1..(end_offset - begin_offset) {
            let cur_op = operations[begin_offset + i];
            if cur_op == last_op {
                last_op_length += 1;
            } else {
                cstr.push_str(&format!("{}", last_op_length));
                cstr.push(last_op as u8 as char);
                last_op = cur_op;
                last_op_length = 1;
            }
        }
        cstr.push_str(&format!("{}", last_op_length));
        cstr.push(last_op as u8 as char);
        cstr
    }

    #[allow(dead_code)]
    fn matching(
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
        let operations = cigar.operations;

        let begin_offset = cigar.begin_offset as isize + offset as isize;
        let end_offset = cigar.end_offset as isize - offset as isize;

        for i in begin_offset..end_offset {
            match operations[i as usize] as u8 as char {
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

#[cfg(test)]
mod tests {
    use super::*;
    use WfaOp::*;

    const PATTERN: &[u8] = b"AGCTAGTGTCAATGGCTACTTTTCAGGTCCT";
    const TEXT: &[u8] = b"AACTAAGTGTCGGTGGCTACTATATATCAGGTCCT";

    fn raw_cigar_string(aligner: &WFAligner) -> String {
        String::from_utf8(aligner.cigar_operations()).unwrap()
    }

    #[test]
    fn test_reap_preserves_aligner_for_reuse() {
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .edit()
            .build();

        aligner.reap();

        let first_result = aligner.align_end_to_end(PATTERN, TEXT);
        let first_cigar = aligner.cigar_operations();
        assert_eq!(first_result.status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(first_result.score, 7);
        assert_eq!(first_cigar, b"MXMMMIMMMMMXXMMMMMMMMIMIMIMMMMMMMMM");

        aligner.reap();
        aligner.reap();

        let second_pattern = b"TCTTTACTCGCGCGTTGGAGAAATACAATAGT";
        let second_text = b"TCTATACTGCGCGTTTGGAGAAATAAAATAGT";
        let second_result = aligner.align_end_to_end(second_pattern, second_text);

        assert_eq!(second_result.status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(second_result.score, 4);
        assert_eq!(
            aligner.cigar_operations(),
            b"MMMXMMMMDMMMMMMMIMMMMMMMMMXMMMMMM"
        );
    }

    #[test]
    fn test_aligner_indel() {
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .indel()
            .build();
        let result = aligner.align_end_to_end(PATTERN, TEXT);
        assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(result.score, 10);
        assert!(!result.dropped);
        assert!(result.null_steps >= 0);
        assert!(result.memory_used > 0);
        assert_eq!(aligner.score(), 10);
        assert_eq!(aligner.cigar_string(None), "1M1I1D3M1I5M2I2D8M1I1M1I1M1I9M");
        let (a, b, c) = aligner.matching(PATTERN, TEXT, None);
        assert_eq!(
            format!("{}\n{}\n{}", a, b, c),
            "A-GCTA-GTGTC--AATGGCTACT-T-T-TCAGGTCCT\n|  ||| |||||    |||||||| | | |||||||||\nAA-CTAAGTGTCGG--TGGCTACTATATATCAGGTCCT"
        );
    }

    #[test]
    fn test_aligner_edit() {
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .edit()
            .build();
        let result = aligner.align_end_to_end(PATTERN, TEXT);
        assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), 7);
        assert_eq!(aligner.cigar_string(None), "1M1X3M1I5M2X8M1I1M1I1M1I9M");
        let (a, b, c) = aligner.matching(PATTERN, TEXT, None);
        assert_eq!(
            format!("{}\n{}\n{}", a, b, c),
            "AGCTA-GTGTCAATGGCTACT-T-T-TCAGGTCCT\n| ||| |||||  |||||||| | | |||||||||\nAACTAAGTGTCGGTGGCTACTATATATCAGGTCCT"
        );
    }

    #[test]
    fn test_aligner_gap_linear() {
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .linear(6, 2)
            .build();
        let result = aligner.align_end_to_end(PATTERN, TEXT);
        assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), -20);
        assert_eq!(aligner.cigar_string(None), "1M1I1D3M1I5M2I2D8M1I1M1I1M1I9M");
        let (a, b, c) = aligner.matching(PATTERN, TEXT, None);
        assert_eq!(
            format!("{}\n{}\n{}", a, b, c),
            "A-GCTA-GTGTC--AATGGCTACT-T-T-TCAGGTCCT\n|  ||| |||||    |||||||| | | |||||||||\nAA-CTAAGTGTCGG--TGGCTACTATATATCAGGTCCT"
        );
    }

    #[test]
    fn test_aligner_gap_affine() {
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
            .affine(6, 4, 2)
            .build();
        let result = aligner.align_end_to_end(PATTERN, TEXT);
        assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), -40);
        assert_eq!(aligner.cigar_string(None), "1M1X3M1I5M2X8M3I1M1X9M");
        let (a, b, c) = aligner.matching(PATTERN, TEXT, None);
        assert_eq!(
            format!("{}\n{}\n{}", a, b, c),
            "AGCTA-GTGTCAATGGCTACT---TTTCAGGTCCT\n| ||| |||||  ||||||||   | |||||||||\nAACTAAGTGTCGGTGGCTACTATATATCAGGTCCT"
        );
    }

    #[test]
    fn test_readme_end_to_end() {
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
            .affine(6, 4, 2)
            .build();

        let pattern = b"TCTTTACTCGCGCGTTGGAGAAATACAATAGT";
        let text = b"TCTATACTGCGCGTTTGGAGAAATAAAATAGT";
        let result = aligner.align_end_to_end(pattern, text);
        assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), -24);
        assert_eq!(
            aligner.cigar_operations(),
            b"MMMXMMMMDMMMMMMMIMMMMMMMMMXMMMMMM"
        );
    }

    #[test]
    fn test_affine_with_match_long_gap() {
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
            .affine_with_match(-1, 2, 2, 1)
            .build();

        let pattern = b"ATAATA";
        let text = b"ATACATAAAATA";
        let result = aligner.align_end_to_end(pattern, text);
        assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), -2);
        assert_eq!(raw_cigar_string(&aligner), "MMMIIIIIIMMM");
    }

    #[test]
    fn test_aligner_score_only() {
        let mut aligner = WFAligner::builder(AlignmentScope::Score, MemoryModel::MemoryLow)
            .affine(6, 4, 2)
            .build();
        let result = aligner.align_end_to_end(PATTERN, TEXT);
        assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), -40);
        assert_eq!(aligner.cigar_string(None), "");
        let (a, b, c) = aligner.matching(PATTERN, TEXT, None);
        assert_eq!(format!("{}\n{}\n{}", a, b, c), "\n\n");
    }

    #[test]
    fn test_aligner_gap_affine_2pieces() {
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .affine2p(6, 2, 2, 4, 1)
            .build();
        let result = aligner.align_end_to_end(PATTERN, TEXT);
        assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), -34);
        assert_eq!(aligner.cigar_string(None), "1M1X3M1I5M2X8M1I1M1I1M1I9M");
        let (a, b, c) = aligner.matching(PATTERN, TEXT, None);
        assert_eq!(
            format!("{}\n{}\n{}", a, b, c),
            "AGCTA-GTGTCAATGGCTACT-T-T-TCAGGTCCT\n| ||| |||||  |||||||| | | |||||||||\nAACTAAGTGTCGGTGGCTACTATATATCAGGTCCT"
        );
    }

    #[test]
    fn test_affine2p_with_match_long_gap() {
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
            .affine2p_with_match(-1, 3, 3, 3, 10, 0)
            .build();

        let pattern = b"TCTATAATAGT";
        let text = b"TCTATACTGCGCGTTTGGAGAAATAAAATAGT";
        let result = aligner.align_end_to_end(pattern, text);
        assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), 1);
        assert_eq!(aligner.cigar_string(None), "6M21I5M");
    }

    #[test]
    fn test_affine2p_with_zero_open() {
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
            .affine2p_with_match(-1, 3, 0, 4, 0, 10)
            .build();

        let pattern = b"TCTATAATAGT";
        let text = b"TCTATACTGCGCGTTTGGAGAAATAAAATAGT";
        let result = aligner.align_end_to_end(pattern, text);
        assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), -73);
        assert_eq!(
            raw_cigar_string(&aligner),
            "MMMMMMIIIIIIIIIIIIMIIIIMMIIIIIMM"
        );
    }

    #[test]
    fn test_linear_and_affine_zero_open_score_equivalence() {
        let pattern = b"ATAATA";
        let text = b"ATACATAAAATA";

        let mut affine_aligner =
            WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
                .affine_with_match(-1, 2, 0, 1)
                .build();
        let result = affine_aligner.align_end_to_end(pattern, text);
        assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(affine_aligner.score(), 0);

        let mut linear_aligner =
            WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
                .linear_with_match(-1, 2, 1)
                .build();
        let result = linear_aligner.align_end_to_end(pattern, text);
        assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(linear_aligner.score(), 0);
    }

    #[test]
    fn test_aligner_span_1() {
        let pattern = b"AATTTAAGTCTAGGCTACTTTC";
        let text = b"CCGACTACTACGAAATTTAAGTATAGGCTACTTTCCGTACGTACGTACGT";
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .affine2p(8, 4, 2, 24, 1)
            .build();
        let result = aligner.align_ends_free(pattern, 0, 0, text, 0, text.len() as i32);
        assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
        let ((xstart, xend), (ystart, yend)) = aligner.get_alignment_span();
        assert_eq!(ystart, 13);
        assert_eq!(yend, 35);
        assert_eq!(xstart, 0);
        assert_eq!(xend, 22);
    }

    #[test]
    fn test_aligner_span_2() {
        let pattern = b"GGGATCCCCGAAAAAGCGGGTTTGGCAAAAGCAAATTTCCCGAGTAAGCAGGCAGAGATCGCGCCAGACGCTCCCCAGAGCAGGGCGTCATGCACAAGAAAGCTTTGCACTTTGCGAACCAACGATAGGTGGGGGTGCGTGGAGGATGGAACACGGACGGCCCGGCTTGCTGCCTTCCCAGGCCTGCAGTTTGCCCATCCACGTCAGGGCCTCAGCCTGGCCGAAAGAAAGAAATGGTCTGTGATCCCCC";
        let text = b"AGCAGGGCGTCATGCACAAGAAAGCTTTGCACTTTGCGAACCAACGATAGGTGGGGGTGCGTGGAGGATGGAACACGGACGGCCCGGCTTGCTGCCTTCCCAGGCCTGCAGTTTGCCCATCCACGTCAGGGCCTCAGCCTGGCCGAAAGAAAGAAATGGTCTGTGATCCCCCCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCATTCCCGGCTACAAGGACCCTTCGAGCCCCGTTCGCCGGCCGCGGACCCGGCCCCTCCCTCCCCGGCCGCTAGGGGGCGGGCCCGGATCACAGGACTGGAGCTGGGCGGAGACCCACGCTCGGAGCGGTTGTGAACTGGCAGGCGGTGGGCGCGGCTTCTGTGCCGTGCCCCGGGCACTCAGTCTTCCAACGGGGCCCCGGAGTCGAAGACAGTTCTAGGGTTCAGGGAGCGCGGGCGGCTCCTGGGCGGCGCCAGACTGCGGTGAGTTGGCCGGCGTGGGCCACCAACCCAATGCAGCCCAGGGCGGCGGCACGAGACAGAACAACGGCGAACAGGAGCAGGGAAAGCGCCTCCGATAGGCCAGGCCTAGGGACCTGCGGGGAGAGGGCGAGGTCAACACCCGGCATGGGCCTCTGATTGGCTCCTGGGACTCGCCCCGCCTACGCCCATAGGTGGGCCCGCACTCTTCCCTGCGCCCCGCCCCCGCCCCAACAGCCT";
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .affine2p(8, 4, 2, 24, 1)
            .with_heuristics(Heuristics::none())
            .build();
        let result = aligner.align_ends_free(pattern, 0, 0, text, 0, text.len() as i32);
        assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
        let ((xstart, xend), (ystart, yend)) = aligner.get_alignment_span();

        assert_eq!(ystart, 0);
        assert_eq!(yend, 172);
        assert_eq!(xstart, 78);
        assert_eq!(xend, 250);
    }

    #[test]
    fn test_aligner_ends_free_global() {
        let pattern = b"AATTTAAGTCTAGGCTACTTTC";
        let text = b"CCGACTACTACGAAATTTAAGTATAGGCTACTTTCCGTACGTACGTACGT";
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .affine(6, 4, 2)
            .build();
        let result = aligner.align_ends_free(pattern, 0, 0, text, 0, text.len() as i32);
        assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), -36);
        assert_eq!(aligner.cigar_string(None), "13I9M1X12M15I");
        let (a, b, c) = aligner.matching(pattern, text, None);
        assert_eq!(
            format!("{}\n{}\n{}", a, b, c),
            "-------------AATTTAAGTCTAGGCTACTTTC---------------\n             ||||||||| ||||||||||||               \nCCGACTACTACGAAATTTAAGTATAGGCTACTTTCCGTACGTACGTACGT"
        );
    }

    #[test]
    fn test_ends_free_with_match_penalties() {
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
            .affine_with_match(-1, 3, 2, 1)
            .build();

        let pattern = b"CGCGTTTGGAGAA";
        let text = b"TCTATACTGCGCGTTTGGAGAAATAAAATAGT";
        let pattern_size = pattern.len() as i32;
        let text_size = text.len() as i32;
        let result = aligner.align_ends_free(
            pattern,
            pattern_size,
            pattern_size,
            text,
            text_size,
            text_size,
        );
        assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), 13);
        assert_eq!(aligner.cigar_string(None), "9I13M10I");

        let pattern = b"TCTATACTGCGCGTTTGGAGAAATAAAATAGT";
        let text = b"CGCGTTTGGAGAA";
        let pattern_size = pattern.len() as i32;
        let text_size = text.len() as i32;
        let result = aligner.align_ends_free(
            pattern,
            pattern_size,
            pattern_size,
            text,
            text_size,
            text_size,
        );
        assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), 13);
        assert_eq!(aligner.cigar_string(None), "9D13M10D");
    }

    #[test]
    fn test_ends_free_shift() {
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
            .affine_with_match(-1, 3, 2, 1)
            .build();

        let pattern = b"TATATTTTTTTTGGAGAAATAAAATA";
        let text = b"TCTATATTTTTTTTTGGAGAAATAAAATAGT";
        let pattern_size = pattern.len() as i32;
        let text_size = text.len() as i32;
        let result = aligner.align_ends_free(
            pattern,
            pattern_size,
            pattern_size,
            text,
            text_size,
            text_size,
        );
        assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(
            raw_cigar_string(&aligner),
            "IIMMMMMMMMMMMMIMMMMMMMMMMMMMMII"
        );
    }

    #[test]
    fn test_aligner_ends_free_right_extent() {
        let pattern = b"AATTTAAGTCTGCTACTTTCACGCAGCT";
        let text = b"AATTTCAGTCTGGCTACTTTCACGTACGATGACAGACTCT";
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .affine(6, 4, 2)
            .build();
        let result =
            aligner.align_ends_free(pattern, 0, pattern.len() as i32, text, 0, text.len() as i32);
        assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), -24);
        assert_eq!(aligner.cigar_string(None), "5M1X6M1I11M4D1M15I");
        let (a, b, c) = aligner.matching(pattern, text, None);
        assert_eq!(
            format!("{}\n{}\n{}", a, b, c),
            "AATTTAAGTCTG-CTACTTTCACGCAGCT---------------\n||||| |||||| |||||||||||    |               \nAATTTCAGTCTGGCTACTTTCACG----TACGATGACAGACTCT"
        );
    }

    #[test]
    fn test_aligner_extension_trims_to_maximal_scoring_prefix() {
        let pattern = b"AATTTAAGTCTGCTACTTTCACGCAGCT";
        let text = b"AATTTCAGTCTGGCTACTTTCACGTACGATGACAGACTCT";

        let mut ends_free_aligner =
            WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
                .affine(6, 4, 2)
                .build();
        let ends_free_result = ends_free_aligner.align_ends_free(
            pattern,
            0,
            pattern.len() as i32,
            text,
            0,
            text.len() as i32,
        );
        assert_eq!(ends_free_result.status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(ends_free_aligner.score(), -24);
        assert_eq!(ends_free_aligner.cigar_string(None), "5M1X6M1I11M4D1M15I");

        let mut extension_aligner =
            WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
                .affine(6, 4, 2)
                .build();
        let extension_result = extension_aligner.align_extension(pattern, text);
        assert_eq!(extension_result.status, AlignmentStatus::StatusAlgPartial);
        assert_eq!(extension_aligner.score(), 10);
        assert_eq!(extension_aligner.cigar_string(None), "5M1X6M1I11M");
        assert_eq!(extension_aligner.cigar_score(), -12);

        let alignment = extension_aligner.get_alignment();
        assert_eq!(alignment.score, 10);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.xend, 23);
        assert_eq!(alignment.ystart, 0);
        assert_eq!(alignment.yend, 24);

        let ((xstart, xend), (ystart, yend)) = extension_aligner.get_alignment_span();
        assert_eq!((xstart, xend), (0, 23));
        assert_eq!((ystart, yend), (0, 24));
    }

    #[test]
    fn test_aligner_extension_empty_prefix_has_zero_span() {
        // No positive-scoring extension exists, so WFA2 trims the entire alignment away. The
        // CIGAR ends up empty, and the span must stay consistent with that (an empty `(0, 0)`
        // span) rather than reflecting the stale wavefront end position.
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .affine(6, 4, 2)
            .build();

        let result = aligner.align_extension(b"AAAAAAAA", b"TTTTTTTT");
        assert_eq!(result.status, AlignmentStatus::StatusAlgPartial);
        assert_eq!(aligner.cigar_string(None), "");

        let alignment = aligner.get_alignment();
        assert!(alignment.operations.is_empty());
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.xend, 0);
        assert_eq!(alignment.ystart, 0);
        assert_eq!(alignment.yend, 0);

        assert_eq!(aligner.get_alignment_span(), ((0, 0), (0, 0)));
    }

    #[test]
    fn test_extension_alignment_span_from_ops() {
        // Anchored at the origin: leading and trailing indels extend the span (unlike the
        // ends-free/local span, which strips them).
        assert_eq!(extension_alignment_span_from_ops(b""), ((0, 0), (0, 0)));
        assert_eq!(extension_alignment_span_from_ops(b"MMM"), ((0, 3), (0, 3)));
        assert_eq!(extension_alignment_span_from_ops(b"IMMM"), ((0, 3), (0, 4)));
        assert_eq!(extension_alignment_span_from_ops(b"MMMII"), ((0, 3), (0, 5)));
        assert_eq!(extension_alignment_span_from_ops(b"DDMMX"), ((0, 5), (0, 3)));
    }

    #[test]
    fn test_aligner_extension_supports_score_scope() {
        let pattern = b"AATTTAAGTCTGCTACTTTCACGCAGCT";
        let text = b"AATTTCAGTCTGGCTACTTTCACGTACGATGACAGACTCT";
        let mut aligner = WFAligner::builder(AlignmentScope::Score, MemoryModel::MemoryHigh)
            .affine(6, 4, 2)
            .build();

        let result = aligner.align_extension(pattern, text);
        assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), -24);
    }

    #[test]
    #[should_panic(expected = "Extension alignment is not supported with MemoryUltraLow")]
    fn test_aligner_extension_rejects_ultralow_memory() {
        let mut aligner =
            WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryUltraLow)
                .affine(6, 4, 2)
                .build();

        aligner.align_extension(b"ACGT", b"ACGT");
    }

    #[test]
    #[should_panic(expected = "Ends-free alignment is not supported with MemoryUltraLow")]
    fn test_aligner_ends_free_rejects_ultralow_memory() {
        let mut aligner =
            WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryUltraLow)
                .affine(6, 4, 2)
                .build();

        aligner.align_ends_free(b"ACGT", 0, 0, b"ACGT", 0, 1);
    }

    #[test]
    fn test_aligner_ends_free_left_extent() {
        let pattern = b"CTTTCACGTACGTGACAGTCTCT";
        let text = b"AATTTCAGTCTGGCTACTTTCACGTACGATGACAGACTCT";
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .affine(6, 4, 2)
            .build();
        let result = aligner.align_ends_free(pattern, 0, 0, text, 0, 0);
        assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), -48);
        assert_eq!(aligner.cigar_string(None), "16I12M1I6M1X4M");
        let (a, b, c) = aligner.matching(pattern, text, None);
        assert_eq!(
            format!("{}\n{}\n{}", a, b, c),
            "----------------CTTTCACGTACG-TGACAGTCTCT\n                |||||||||||| |||||| ||||\nAATTTCAGTCTGGCTACTTTCACGTACGATGACAGACTCT"
        );

        // Leading 16I is stripped (`ystart` = 16); there are no trailing indels, so the span
        // runs to the end of both sequences.
        assert_eq!(aligner.get_alignment_span(), ((0, 23), (16, 40)));
    }

    #[test]
    fn test_aligner_ends_free_right_overlap() {
        let pattern = b"CGCGTCTGACTGACTGACTAAACTTTCATGTACCTGACA";
        let text = b"AAACTTTCACGTACGTGACATATAGCGATCGATGACT";
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .affine(6, 4, 2)
            .build();
        let result = aligner.align_ends_free(pattern, 0, 0, text, 0, 0);
        assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), -92);
        assert_eq!(aligner.cigar_string(None), "19D9M1X4M1X5M17I");
        let (a, b, c) = aligner.matching(pattern, text, None);
        assert_eq!(
            format!("{}\n{}\n{}", a, b, c),
            "CGCGTCTGACTGACTGACTAAACTTTCATGTACCTGACA-----------------\n                   ||||||||| |||| |||||                 \n-------------------AAACTTTCACGTACGTGACATATAGCGATCGATGACT"
        );

        // The span is symmetric: leading 19D and trailing 17I are both stripped, so it covers
        // only the aligned core. `yend` stops at the last M/X column (20), not the full text.
        assert_eq!(aligner.get_alignment_span(), ((19, 39), (0, 20)));
    }

    #[test]
    fn test_ends_free_span_excludes_trailing_insertions() {
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .linear_with_match(-1, 1, 1)
            .build();

        let pattern = b"A";
        let text = b"ACG";
        let result = aligner.align_ends_free(pattern, 0, 0, text, 0, 2);
        assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), 1);
        assert_eq!(raw_cigar_string(&aligner), "MII");

        // The trailing II remain in the active CIGAR ("MII"), but the span strips trailing
        // indels, so it covers only the single matched column.
        assert_eq!(aligner.get_alignment_span(), ((0, 1), (0, 1)));
    }

    #[test]
    fn test_clipping_score() {
        let text_lf = b"AAGGAGCTGAGAATTGTTCTTCCAGATACCTTTCCGACCTCTTCTTGGTT";
        let text_rf = b"GGAGTGCAGTGGTGCAATCTTGGCTCACTACAACCTCCGCATCCTGGGTT";

        let pattern_lf = b"AAGGAGCTGAGAATTGTTCGTCCAGATACCTTTCCGACCTCTTCTTGGTT";
        let pattern_rf = b"GGAGTGCAGTGGTGCAATCTTGGCTCACTACAACCTCTGCATCCTGGGTT";

        let motif = b"ATTT";

        let text = [text_lf, &motif.repeat(10)[..], text_rf].concat();
        let pattern = [pattern_lf, &motif.repeat(8)[..], pattern_rf].concat();

        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .affine2p(8, 4, 2, 24, 1)
            .build();
        let result = aligner.align_end_to_end(&pattern, &text);

        assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), -36);
        assert_eq!(aligner.cigar_string(None), "19M1X62M8I37M1X12M");
        let (a, b, c) = aligner.matching(&pattern, &text, None);
        assert_eq!(
            format!("{}\n{}\n{}", a, b, c),
            "AAGGAGCTGAGAATTGTTCGTCCAGATACCTTTCCGACCTCTTCTTGGTTATTTATTTATTTATTTATTTATTTATTTATTT--------GGAGTGCAGTGGTGCAATCTTGGCTCACTACAACCTCTGCATCCTGGGTT\n||||||||||||||||||| ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||        ||||||||||||||||||||||||||||||||||||| ||||||||||||\nAAGGAGCTGAGAATTGTTCTTCCAGATACCTTTCCGACCTCTTCTTGGTTATTTATTTATTTATTTATTTATTTATTTATTTATTTATTTGGAGTGCAGTGGTGCAATCTTGGCTCACTACAACCTCCGCATCCTGGGTT"
        );
        assert_eq!(aligner.cigar_score(), -36);
        assert_eq!(aligner.cigar_score_clipped(50), -20);
        assert_eq!(aligner.cigar_string(Some(50)), "32M8I");
        let (a, b, c) = aligner.matching(&pattern, &text, Some(50));
        assert_eq!(
            format!("{}\n{}\n{}", a, b, c),
            "ATTTATTTATTTATTTATTTATTTATTTATTT--------\n||||||||||||||||||||||||||||||||        \nATTTATTTATTTATTTATTTATTTATTTATTTATTTATTT"
        );

        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .indel()
            .with_heuristics(Heuristics::none())
            .build();
        let result = aligner.align_end_to_end(&pattern, &text);
        assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), 12);
        assert_eq!(aligner.cigar_score(), 12);
        assert_eq!(aligner.cigar_score_clipped(19), 10);
        assert_eq!(aligner.cigar_score_clipped(0), 12);
    }

    #[test]
    fn test_memory_modes() {
        let expected_cigar = "1M1X3M1I5M2X8M3I1M1X9M";
        let expected_matching = "AGCTA-GTGTCAATGGCTACT---TTTCAGGTCCT\n| ||| |||||  ||||||||   | |||||||||\nAACTAAGTGTCGGTGGCTACTATATATCAGGTCCT";
        let expected_score = -48;

        struct Test {
            memory_mode: MemoryModel,
        }

        let tests = vec![
            Test {
                memory_mode: MemoryModel::MemoryHigh,
            },
            Test {
                memory_mode: MemoryModel::MemoryMed,
            },
            Test {
                memory_mode: MemoryModel::MemoryLow,
            },
            // Test {
            //     memory_mode: MemoryModel::MemoryUltraLow,
            // },
        ];

        for test in tests {
            let mut aligner = WFAligner::builder(AlignmentScope::Alignment, test.memory_mode)
                .affine2p(8, 4, 2, 24, 1)
                .build();
            let result = aligner.align_end_to_end(PATTERN, TEXT);
            assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
            assert_eq!(aligner.score(), expected_score);
            assert_eq!(aligner.cigar_score(), expected_score);
            assert_eq!(aligner.cigar_score_clipped(0), expected_score);
            assert_eq!(aligner.cigar_string(None), expected_cigar);
            let (a, b, c) = aligner.matching(PATTERN, TEXT, None);
            assert_eq!(format!("{}\n{}\n{}", a, b, c), expected_matching);
        }
    }

    #[test]
    fn test_set_heuristics_replaces_configuration() {
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .affine(6, 4, 2)
            .build();
        let combined = Heuristics::new(3)
            .with_adaptive(AdaptiveHeuristic::WfMash {
                min_wavefront_length: 1,
                max_distance_threshold: 2,
            })
            .with_drop(DropHeuristic::XDrop { xdrop: 10 })
            .with_band(BandHeuristic::Adaptive { min_k: 1, max_k: 2 });
        aligner.set_heuristics(combined);
        assert_eq!(aligner.get_heuristics(), combined);

        let replacement = Heuristics::banded_static(1, 2);
        aligner.set_heuristics(replacement);
        assert_eq!(aligner.get_heuristics(), replacement);
    }

    #[test]
    fn test_resource_limits_builder_and_setters() {
        let initial_limits = ResourceLimits::new(64, 1_048_576, 2_097_152, 1, 64);
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .with_resource_limits(initial_limits)
            .affine(6, 4, 2)
            .build();
        assert_eq!(aligner.get_resource_limits(), initial_limits);

        aligner.set_max_alignment_steps(128);
        aligner.set_max_memory(2_097_152, 4_194_304);
        aligner.set_max_num_threads(2);
        aligner.set_min_offsets_per_thread(32);

        assert_eq!(
            aligner.get_resource_limits(),
            ResourceLimits {
                max_alignment_steps: 128,
                max_memory_resident: 2_097_152,
                max_memory_abort: 4_194_304,
                max_num_threads: 2,
                min_offsets_per_thread: 32,
            }
        );

        let aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .with_max_alignment_steps(256)
            .with_max_memory(4_194_304, 8_388_608)
            .with_max_num_threads(1)
            .with_min_offsets_per_thread(128)
            .affine(6, 4, 2)
            .build();
        assert_eq!(
            aligner.get_resource_limits(),
            ResourceLimits {
                max_alignment_steps: 256,
                max_memory_resident: 4_194_304,
                max_memory_abort: 8_388_608,
                max_num_threads: 1,
                min_offsets_per_thread: 128,
            }
        );
    }

    #[test]
    fn test_max_alignment_steps_limit() {
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .with_max_alignment_steps(1)
            .edit()
            .build();

        let result = aligner.align_end_to_end(PATTERN, TEXT);
        assert_eq!(result.status, AlignmentStatus::StatusMaxStepsReached);
        assert!(!result.dropped);
        assert!(result.null_steps >= 0);
    }

    #[test]
    #[should_panic(expected = "No valid alignment span is available")]
    fn test_get_alignment_span_rejects_missing_cigar_end() {
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .with_max_alignment_steps(1)
            .edit()
            .build();

        let result = aligner.align_end_to_end(PATTERN, TEXT);
        assert_eq!(result.status, AlignmentStatus::StatusMaxStepsReached);
        aligner.get_alignment_span();
    }

    #[test]
    #[should_panic(expected = "No valid alignment span is available")]
    fn test_get_alignment_rejects_missing_cigar_end() {
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .with_max_alignment_steps(1)
            .edit()
            .build();

        let result = aligner.align_end_to_end(PATTERN, TEXT);
        assert_eq!(result.status, AlignmentStatus::StatusMaxStepsReached);
        aligner.get_alignment();
    }

    #[test]
    #[should_panic(expected = "Cannot get alignment when AlignmentScope is Score")]
    fn test_get_alignment_rejects_score_scope() {
        let aligner = WFAligner::builder(AlignmentScope::Score, MemoryModel::MemoryHigh)
            .edit()
            .build();

        aligner.get_alignment();
    }

    #[test]
    #[should_panic(expected = "Cannot get alignment span when AlignmentScope is Score")]
    fn test_get_alignment_span_rejects_score_scope() {
        let aligner = WFAligner::builder(AlignmentScope::Score, MemoryModel::MemoryHigh)
            .edit()
            .build();

        aligner.get_alignment_span();
    }

    #[test]
    fn ultralow_memory_default_heuristic_can_be_unattainable() {
        let read = b"GCTGCTACTGGGGTGTCCCCTCTCAAAGGACAAACCCAGGATCTACAGATGTGTGTGCTAAGCCATGTATGCACATGCACGTGTGTGTGTATATATTTAACCTATCTGTATATATGTATTATGTAAACATGAGTTCCTGCTGGCATATCTGACTATAACTGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACATGAGTTCCCTGCTGGCATATCTGACTATAACTGACCACCTCACAGTCCATTCTGATCTCTATATATGTATTATGTAAACATGAGTTCCTACTGGCATATCTGACTATAACTGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATTATGTAAACATGAGTTCCCTGCTGGCATATCTGATTATAACTGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATTATGTAAACATGAGTTCCTACTGGCATATCTGACTATAACCGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACACGAGTTCCTACTGGCATATCTGACTATAACTGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACACGAGTTCCTGCTGGCATATCTGACTATAACTGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATAATATATATTATATATGGACCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACATGAGTTCCTGCTGGCATATCTGACTATAACTGACCACCTCAGGGTCTATTCTGATCTGTATATATGTATAATATATATTATATATGGACCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACATGAGTTCCTGCTGGCATATCTGATTATAACCGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATTATGTAAACATGAGTTCCTACTGGCATATCTGACTATAACCGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACATGAGTTCCTACTGGCATATCTGACTATAACTGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACATGAGTTCCTACTGGCATATCTGACTATAACTGACCACCTCAGGATCCATTCTGATCTGTATATATGTATAATATATATTATATATGGACCTCAGGGTCCATTCTGATCTGTATATATGTATAATATATATTATATATGGACCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACATGAGTTCCTGGCTGGCATATCTGATTATAACCGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACACGAGTTCCTACTGGCATATCTGACTATAACTGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACACGAGTTCCTGCTGGCATATCTGATTATAACCGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATAATATATATTATATATGGACCTCAGGGTCCCCGCTGGCTTTTCCATGACTTCCTTATCCAGCTGTGAGAACCCTGACTCTTACTACCCATACTGTATTGACTTATTT";
        let allele = b"GCTGCTACTGGGGTGTCCCCTCTCAAAGGACAAACCCAGGATCTACAGATGTGTGTGCTAAGCCATGTATGCACACGCACGTGTGTGTGTATATATTTAACCTATCTGTATATATGTATTATGTAAACATGAGTTCCTGCTGGCATATCTGACTATAACTGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACACGACTTCCTACTGGCATATCTGACTGTAACCGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACATGATTTCCTACTGGCATATCTGACTATAACTGACCACCTCAGGGTTCATTCCGATCTGTATATAAGTATCATGTAAACACGAGTTCCTGCTGGCATATCTGACTGTAACCGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACACGAGTTCCTGCTGGCATATCTGACTATAACCGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACATGAGTTCCTACTGGCATATCTGACTATAACTGACCACCTCAGGGTCCATTCTGATCTGTATGTATGTATCATGTAAACACGAGTTCCTACTGGCATATCTGACTATAACTGACCACCTCAGGGTCCATTCCGATCTGTATATAAGTATCATGTAAACACGAGTTCCTGCTGGCATATCTGACTGTAACCGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACACGAGTTCCTGCTGGCATATCTGACTATAACTGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATAATATATATTATATATGGACCTCAGGGTCCATTCTGATCTGCATATATGTATAATATATATTATATATGGACCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACATGAGTTCCTGCTGGCATATCTGTCTATAACCGACCACCTTAGGGTCCATTCTGATCTGTATATATGTATAATATATATTATATATGGTCCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACATGAGTTCCTGCTGGCATATCTGTCTATAACCGACCACCTTAGGGTCCATTCTGATCTGTATATATGTATAATATATATTATATATGGACCTCAGGGTCCATTCTGATCTGCATATATGTATAATATATATTATATATGGTCCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACATGAGTTCCTGCTGGCATATCTGTCTATAACCGACCACCTTAGGGTCCATTCTGATCTGTATATATGTATAATATATATTATATATGGACCTCAGGGTCCCCGCTGGCTTTTCCATGACTTCCTTATCCAGCTGTGAGAACCCTGACTCTTACTACTGTATTGACTTATTTGTGAAACCT";

        let mut aligner =
            WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryUltraLow)
                .with_heuristics(Heuristics::wfa2_default())
                .affine2p(8, 4, 2, 24, 1)
                .build();

        // These are valid sequence inputs, the unattainable result is specific to WFA2's
        // BiWFA path when `wf_adaptive(1, 10, 50)` (i.e. `wfa2_default()`).
        // For this pair, this heuristic prunes enough state that BiWFA reaches an end before it
        // can find a midpoint breakpoint. The reached score is above WFA2's
        // BiWFA recovery threshold, so WFA2 reports `WF_STATUS_UNATTAINABLE`.
        let result = aligner.align_end_to_end(read, allele);
        assert_eq!(result.status, AlignmentStatus::StatusUnattainable);
        assert_eq!(aligner.score(), i32::MIN);

        // Setting a more permissive heuristic allows BiWFA to find a midpoint
        // breakpoint and recover with its regular fallback path.
        aligner.set_heuristics(Heuristics::wf_adaptive(1, 10, 75));
        let result = aligner.align_end_to_end(read, allele);
        assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), -881);

        aligner.set_heuristics(Heuristics::none());
        let result = aligner.align_end_to_end(read, allele);
        assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), -881);
    }

    #[test]
    fn test_get_penalties() {
        let aligner_edit = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
            .edit()
            .build();
        assert_eq!(aligner_edit.get_penalties(), Penalties::Edit);

        let aligner_indel = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
            .indel()
            .build();
        assert_eq!(aligner_indel.get_penalties(), Penalties::Indel);

        let aligner_linear = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
            .linear(12, 24)
            .build();
        assert_eq!(
            aligner_linear.get_penalties(),
            Penalties::Linear {
                match_: 0,
                mismatch: 12,
                indel: 24
            }
        );

        let aligner_affine = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
            .affine(12, 24, 2)
            .build();
        assert_eq!(
            aligner_affine.get_penalties(),
            Penalties::Affine {
                match_: 0,
                mismatch: 12,
                gap_opening: 24,
                gap_extension: 2
            }
        );

        let aligner_affine2p =
            WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
                .affine2p(12, 24, 2, 48, 1)
                .build();
        assert_eq!(
            aligner_affine2p.get_penalties(),
            Penalties::Affine2p {
                match_: 0,
                mismatch: 12,
                gap_opening1: 24,
                gap_extension1: 2,
                gap_opening2: 48,
                gap_extension2: 1
            }
        );

        let aligner_affine_match =
            WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
                .affine_with_match(-5, 12, 24, 2)
                .build();
        assert_eq!(
            aligner_affine_match.get_penalties(),
            Penalties::Affine {
                match_: -5,
                mismatch: 12,
                gap_opening: 24,
                gap_extension: 2
            }
        );
    }

    #[test]
    fn test_builder_pattern() {
        let aligner_edit = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
            .edit()
            .build();
        assert_eq!(aligner_edit.get_penalties(), Penalties::Edit);

        let aligner_affine = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
            .affine(12, 24, 2)
            .with_heuristics(Heuristics::wf_adaptive(100, 10, 50))
            .build();
        assert_eq!(
            aligner_affine.get_penalties(),
            Penalties::Affine {
                match_: 0,
                mismatch: 12,
                gap_opening: 24,
                gap_extension: 2
            }
        );

        let aligner_affine_heuristic =
            WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
                .affine(12, 24, 2)
                .with_heuristics(Heuristics::wf_adaptive(100, 10, 50))
                .build();
        assert_eq!(
            aligner_affine_heuristic.get_penalties(),
            Penalties::Affine {
                match_: 0,
                mismatch: 12,
                gap_opening: 24,
                gap_extension: 2
            }
        );

        let aligner_affine2p =
            WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
                .affine2p(12, 24, 2, 48, 1)
                .build();
        assert_eq!(
            aligner_affine2p.get_penalties(),
            Penalties::Affine2p {
                match_: 0,
                mismatch: 12,
                gap_opening1: 24,
                gap_extension1: 2,
                gap_opening2: 48,
                gap_extension2: 1
            }
        );

        let aligner_linear = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
            .linear_with_match(-5, 12, 24)
            .build();
        assert_eq!(
            aligner_linear.get_penalties(),
            Penalties::Linear {
                match_: -5,
                mismatch: 12,
                indel: 24
            }
        );
    }

    #[test]
    fn test_get_and_decode_sam_cigar() {
        let pattern = b"TCTTTACTCTT";
        let text = b"TCTTTACTCTT";
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .affine(4, 6, 2)
            .build();

        let result = aligner.align_end_to_end(pattern, text);
        assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);

        let sam_cigar_buffer = aligner.get_sam_cigar(true);
        assert!(
            !sam_cigar_buffer.is_empty(),
            "SAM CIGAR buffer should not be empty"
        );

        let decoded_cigar = WFAligner::decode_sam_cigar(&sam_cigar_buffer);

        // Expected result for identical sequences (11 matches), The raw buffer encodes length << 4 | op_code. '=' is op_code 7. So, 11= should be encoded as (11 << 4) | 7 = 176 | 7 = 183
        let expected_raw_buffer = vec![183]; // 11=
        assert_eq!(
            sam_cigar_buffer, expected_raw_buffer,
            "Raw SAM CIGAR buffer mismatch"
        );

        let expected_decoded_cigar = vec![(11, '=')]; // 11 matches ('=' because show_mismatches=true)
        assert_eq!(
            decoded_cigar, expected_decoded_cigar,
            "Decoded SAM CIGAR mismatch"
        );

        // Test with show_mismatches = false
        let sam_cigar_buffer_m = aligner.get_sam_cigar(false);
        // 'M' is op_code 0. (11 << 4) | 0 = 176
        let expected_raw_buffer_m = vec![176]; // 11M
        assert_eq!(
            sam_cigar_buffer_m, expected_raw_buffer_m,
            "Raw SAM CIGAR buffer mismatch (M)"
        );

        let decoded_cigar_m = WFAligner::decode_sam_cigar(&sam_cigar_buffer_m);
        let expected_decoded_cigar_m: Vec<CigarOp> = vec![(11, 'M')]; // 11 matches ('M')
        assert_eq!(
            decoded_cigar_m, expected_decoded_cigar_m,
            "Decoded SAM CIGAR mismatch (M)"
        );

        let pattern_diff = b"TCTTTACTCTT";
        let text_diff = b"TCTTTACTATT";
        let mut aligner_diff =
            WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
                .affine(4, 6, 2)
                .build();
        let result_diff = aligner_diff.align_end_to_end(pattern_diff, text_diff);
        assert_eq!(result_diff.status, AlignmentStatus::StatusAlgCompleted);

        let sam_cigar_buffer_diff = aligner_diff.get_sam_cigar(true);

        let expected_raw_diff = vec![135, 24, 39];
        assert_eq!(
            sam_cigar_buffer_diff, expected_raw_diff,
            "Raw SAM CIGAR buffer mismatch (diff)"
        );

        let decoded_cigar_diff = WFAligner::decode_sam_cigar(&sam_cigar_buffer_diff);
        let expected_decoded_diff: Vec<CigarOp> = vec![(8, '='), (1, 'X'), (2, '=')];
        assert_eq!(
            decoded_cigar_diff, expected_decoded_diff,
            "Decoded SAM CIGAR mismatch (diff)"
        );

        // Test with show_mismatches = false
        let sam_cigar_buffer_diff_m = aligner_diff.get_sam_cigar(false);
        // Expected: 11M => (11<<4)|0 = 176
        let expected_raw_diff_m = vec![176];
        assert_eq!(
            sam_cigar_buffer_diff_m, expected_raw_diff_m,
            "Raw SAM CIGAR buffer mismatch (diff, M)"
        );

        let decoded_cigar_diff_m = WFAligner::decode_sam_cigar(&sam_cigar_buffer_diff_m);
        let expected_decoded_diff_m: Vec<CigarOp> = vec![(11, 'M')];
        assert_eq!(
            decoded_cigar_diff_m, expected_decoded_diff_m,
            "Decoded SAM CIGAR mismatch (diff, M)"
        );
    }

    #[test]
    fn test_get_heuristics_round_trips_combined_categories() {
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .affine(12, 24, 2)
            .build();
        assert_eq!(aligner.get_heuristics(), Heuristics::none());

        let empty_with_custom_steps = Heuristics::new(10);
        aligner.set_heuristics(empty_with_custom_steps);
        assert_eq!(aligner.get_heuristics(), empty_with_custom_steps);

        let combined = Heuristics::new(5)
            .with_adaptive(AdaptiveHeuristic::WfAdaptive {
                min_wavefront_length: 5,
                max_distance_threshold: 25,
            })
            .with_drop(DropHeuristic::XDrop { xdrop: 15 })
            .with_band(BandHeuristic::Static {
                min_k: 5,
                max_k: 20,
            });
        aligner.set_heuristics(combined);
        assert_eq!(aligner.get_heuristics(), combined);

        let aligner_with_heuristics =
            WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
                .affine(12, 24, 2)
                .with_heuristics(Heuristics::wf_adaptive(100, 10, 50))
                .build();
        assert_eq!(
            aligner_with_heuristics.get_heuristics(),
            Heuristics::wf_adaptive(100, 10, 50)
        );
        assert_eq!(
            aligner_with_heuristics.get_penalties(),
            Penalties::Affine {
                match_: 0,
                mismatch: 12,
                gap_opening: 24,
                gap_extension: 2
            }
        );
    }

    #[test]
    fn test_heuristics_none_clears_configuration() {
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .affine(12, 24, 2)
            .with_heuristics(
                Heuristics::new(5)
                    .with_adaptive(AdaptiveHeuristic::WfAdaptive {
                        min_wavefront_length: 5,
                        max_distance_threshold: 25,
                    })
                    .with_drop(DropHeuristic::ZDrop { zdrop: 15 })
                    .with_band(BandHeuristic::Adaptive {
                        min_k: 5,
                        max_k: 20,
                    }),
            )
            .build();

        aligner.set_heuristics(Heuristics::none());
        assert!(aligner.get_heuristics().is_none());
        assert_eq!(aligner.get_heuristics(), Heuristics::none());
    }

    #[test]
    fn test_drop_heuristics_reject_edit_and_indel() {
        assert!(std::panic::catch_unwind(|| {
            let _ = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
                .edit()
                .with_heuristics(Heuristics::xdrop(1, 10))
                .build();
        })
        .is_err());

        assert!(std::panic::catch_unwind(|| {
            let _ = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
                .indel()
                .with_heuristics(Heuristics::zdrop(1, 10))
                .build();
        })
        .is_err());

        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .edit()
            .build();
        assert!(std::panic::catch_unwind(move || {
            aligner.set_heuristics(Heuristics::xdrop(1, 10));
        })
        .is_err());
    }

    #[test]
    fn test_combined_heuristics_alignment_completes() {
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .affine(6, 4, 2)
            .with_heuristics(
                Heuristics::new(1)
                    .with_adaptive(AdaptiveHeuristic::WfAdaptive {
                        min_wavefront_length: 1,
                        max_distance_threshold: 100,
                    })
                    .with_drop(DropHeuristic::XDrop { xdrop: 1_000 })
                    .with_band(BandHeuristic::Static {
                        min_k: -100,
                        max_k: 100,
                    }),
            )
            .build();

        let result = aligner.align_end_to_end(PATTERN, TEXT);
        assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
    }

    #[test]
    fn test_dropped_alignment_trims_to_maximal_scoring_prefix() {
        let pattern = b"AAAAAAAAAACCCCCCCCCCAAAAAAAAAA";
        let text = b"AAAAAAAAAAGGGGGGGGGGAAAAAAAAAA";
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .affine_with_match(-1, 4, 6, 2)
            .with_heuristics(Heuristics::zdrop(1, 0))
            .build();

        let result = aligner.align_end_to_end(pattern, text);
        assert_eq!(result.status, AlignmentStatus::StatusAlgPartial);
        assert!(result.dropped);
        assert_eq!(aligner.score(), 10);
        assert_eq!(aligner.cigar_string(None), "10M");
        assert_eq!(aligner.cigar_score(), 10);
        assert_eq!(aligner.get_alignment_span(), ((0, 10), (0, 10)));
    }

    #[test]
    fn test_alignment_span_from_ops() {
        // Mixed: leading insertions offset the text start, trailing indels do not extend the span.
        assert_eq!(alignment_span_from_ops(b"IIIMMMDDXII"), ((0, 6), (3, 7)));
        // No aligned columns at all -> empty span on both axes.
        assert_eq!(alignment_span_from_ops(b""), ((0, 0), (0, 0)));
        assert_eq!(alignment_span_from_ops(b"DDII"), ((0, 0), (0, 0)));
        assert_eq!(alignment_span_from_ops(b"III"), ((0, 0), (0, 0)));
        assert_eq!(alignment_span_from_ops(b"DDD"), ((0, 0), (0, 0)));
        // Single aligned column.
        assert_eq!(alignment_span_from_ops(b"M"), ((0, 1), (0, 1)));
        // Substitutions advance both pattern and text just like matches.
        assert_eq!(alignment_span_from_ops(b"XXX"), ((0, 3), (0, 3)));
        // Leading deletions offset the pattern start and leading insertions offset the text start.
        assert_eq!(alignment_span_from_ops(b"DDDMM"), ((3, 5), (0, 2)));
        assert_eq!(alignment_span_from_ops(b"IIIMM"), ((0, 2), (3, 5)));
        // Trailing indels after the last aligned column do not extend the span.
        assert_eq!(alignment_span_from_ops(b"MMIID"), ((0, 2), (0, 2)));
        // Internal gaps diverge the pattern and text spans.
        assert_eq!(alignment_span_from_ops(b"MMDDMM"), ((0, 6), (0, 4)));
        assert_eq!(alignment_span_from_ops(b"MMIIMM"), ((0, 4), (0, 6)));
    }

    #[test]
    fn test_get_alignment_global() {
        let pattern = b"AGCTAGTGTCAATGGCTACTTTTCAGGTCCT";
        let text = b"AACTAAGTGTCGGTGGCTACTATATATCAGGTCCT";
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .affine(1, 5, 1)
            .build();
        let result = aligner.align_end_to_end(pattern, text);
        assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);

        let alignment = aligner.get_alignment();
        assert_eq!(aligner.score(), -18);

        let expected_ops = vec![
            Match, Subst, Match, Match, Match, Ins, Match, Match, Match, Match, Match, Subst,
            Subst, Match, Match, Match, Match, Match, Match, Match, Match, Ins, Ins, Ins, Match,
            Subst, Match, Match, Match, Match, Match, Match, Match, Match, Match,
        ];

        assert_eq!(alignment.score, -18);
        assert_eq!(alignment.xlen, pattern.len());
        assert_eq!(alignment.ylen, text.len());
        assert_eq!(alignment.operations, expected_ops);

        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.xend, 31);
        assert_eq!(alignment.ystart, 0);
        assert_eq!(alignment.yend, 35);

        let ((xstart, xend), (ystart, yend)) = aligner.get_alignment_span();
        assert_eq!(alignment.xstart, xstart);
        assert_eq!(alignment.xend, xend);
        assert_eq!(alignment.ystart, ystart);
        assert_eq!(alignment.yend, yend);
    }

    #[test]
    fn test_get_alignment_biwfa_global() {
        let pattern = b"AGCTAGTGTCAATGGCTACTTTTCAGGTCCT";
        let text = b"AACTAAGTGTCGGTGGCTACTATATATCAGGTCCT";
        let mut aligner =
            WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryUltraLow)
                .affine(1, 5, 1)
                .build();
        let result = aligner.align_end_to_end(pattern, text);
        assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);

        let alignment = aligner.get_alignment();
        assert_eq!(aligner.score(), -18);
        assert_eq!(aligner.cigar_score(), -18);

        let expected_ops = vec![
            Match, Subst, Match, Match, Match, Ins, Match, Match, Match, Match, Match, Subst,
            Subst, Match, Match, Match, Match, Match, Match, Match, Match, Ins, Ins, Ins, Match,
            Subst, Match, Match, Match, Match, Match, Match, Match, Match, Match,
        ];

        assert_eq!(alignment.score, -18);
        assert_eq!(alignment.xlen, pattern.len());
        assert_eq!(alignment.ylen, text.len());
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.xend, pattern.len());
        assert_eq!(alignment.ystart, 0);
        assert_eq!(alignment.yend, text.len());
        assert_eq!(alignment.operations, expected_ops);

        let ((xstart, xend), (ystart, yend)) = aligner.get_alignment_span();
        assert_eq!(alignment.xstart, xstart);
        assert_eq!(alignment.xend, xend);
        assert_eq!(alignment.ystart, ystart);
        assert_eq!(alignment.yend, yend);
    }

    #[test]
    fn test_get_alignment_biwfa_global_long_recursion() {
        // Sequences long and divergent enough to push BiWFA past its fallback thresholds
        // (MIN_LENGTH = 100, MIN_SCORE = 250), forcing multiple recursive splits. Each split
        // rewrites the C aligner's `wf_forward` sequence bounds, so the reported sequence
        // lengths must come from the values captured at `align` time, not the C struct.
        let bases = [b'A', b'C', b'G', b'T'];
        let mut pattern = Vec::new();
        let mut text = Vec::new();
        let mut state: u64 = 0x9E3779B97F4A7C15;
        let mut next_base = || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            bases[((state >> 33) % 4) as usize]
        };
        for _ in 0..400 {
            pattern.push(next_base());
            text.push(next_base());
        }

        let mut aligner =
            WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryUltraLow)
                .affine(1, 5, 1)
                .build();
        let result = aligner.align_end_to_end(&pattern, &text);
        assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);

        let alignment = aligner.get_alignment();
        assert_eq!(alignment.xlen, pattern.len());
        assert_eq!(alignment.ylen, text.len());
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.xend, pattern.len());
        assert_eq!(alignment.ystart, 0);
        assert_eq!(alignment.yend, text.len());

        let ((xstart, xend), (ystart, yend)) = aligner.get_alignment_span();
        assert_eq!((xstart, xend), (0, pattern.len()));
        assert_eq!((ystart, yend), (0, text.len()));
    }

    #[test]
    fn test_get_alignment_ends_free() {
        let pattern = b"AGTGTCAATGGCTAC";
        let text = b"GGGGGGGGGGAGTGTCAATGGCTACGGGGGGGGGG";
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .affine(1, 5, 1)
            .build();
        let result =
            aligner.align_ends_free(pattern, 0, 0, text, text.len() as i32, text.len() as i32);
        assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
        let alignment = aligner.get_alignment();
        assert_eq!(aligner.score(), 0);

        let expected_ops = vec![
            Ins, Ins, Ins, Ins, Ins, Ins, Ins, Ins, Ins, Ins, Match, Match, Match, Match, Match,
            Match, Match, Match, Match, Match, Match, Match, Match, Match, Match, Ins, Ins, Ins,
            Ins, Ins, Ins, Ins, Ins, Ins, Ins,
        ];

        assert_eq!(alignment.score, 0);
        assert_eq!(alignment.xlen, pattern.len());
        assert_eq!(alignment.ylen, text.len());
        assert_eq!(alignment.operations, expected_ops);

        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.xend, pattern.len());
        assert_eq!(alignment.ystart, 10);
        assert_eq!(alignment.yend, 25);

        let ((xstart, xend), (ystart, yend)) = aligner.get_alignment_span();
        assert_eq!(alignment.xstart, xstart);
        assert_eq!(alignment.xend, xend);
        assert_eq!(alignment.ystart, ystart);
        assert_eq!(alignment.yend, yend);
    }
}
