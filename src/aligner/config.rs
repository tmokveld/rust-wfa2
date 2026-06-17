use crate::wfa2;
use std::fmt;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum MemoryModel {
    MemoryHigh,
    MemoryMed,
    MemoryLow,
    MemoryUltraLow,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum WfaOp {
    Match,
    Subst,
    Ins,
    Del,
}

impl WfaOp {
    pub(crate) fn from_u8(op_char: u8) -> Self {
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
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
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
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
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
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
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
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum DropHeuristic {
    XDrop { xdrop: i32 },
    ZDrop { zdrop: i32 },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum BandHeuristic {
    Static { min_k: i32, max_k: i32 },
    Adaptive { min_k: i32, max_k: i32 },
}

/// WFA2 heuristic configuration.
///
/// Rust intentionally defaults to [`Heuristics::none`], unlike WFA2's C
/// `wavefront_aligner_attr_default`, which enables WF-adaptive pruning with
/// `min_wavefront_length = 10`, `max_distance_threshold = 50`, and
/// `steps_between_cutoffs = 1`. This wrapper keeps exact alignment as the
/// default and requires callers to opt into heuristic pruning explicitly.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
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
    /// Disable WFA2 heuristic pruning.
    ///
    /// This is the Rust wrapper default. It intentionally differs from WFA2's C
    /// default attributes, which enable WF-adaptive pruning.
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

    /// Return the heuristic configuration from WFA2's C default attributes.
    ///
    /// Use this when porting C code that relied on
    /// `wavefront_aligner_attr_default.heuristic`.
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

    pub(crate) fn validate(&self) {
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
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
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

/// Options for recording WFA2's native wavefront plot dump.
///
/// Plotting is intended for debugging and external tooling. It records
/// WFA2's upstream `.plot` text format during alignment; it does not render an
/// image.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PlotOptions {
    /// Total heatmap resolution points used by WFA2's plot recorder.
    ///
    /// Default: `2000`.
    pub resolution_points: i32,
    /// BiWFA recursion level to record. `-1` records the final/subsidiary
    /// alignment, and non-negative values record that recursion level.
    ///
    /// Default: `0`.
    pub align_level: i32,
}

impl Default for PlotOptions {
    fn default() -> Self {
        Self {
            resolution_points: 2000,
            align_level: 0,
        }
    }
}

impl PlotOptions {
    pub fn new(resolution_points: i32, align_level: i32) -> Self {
        validate_plot_options(resolution_points, align_level);
        Self {
            resolution_points,
            align_level,
        }
    }

    pub fn final_alignment() -> Self {
        Self {
            align_level: -1,
            ..Self::default()
        }
    }

    pub fn at_recursion_level(level: i32) -> Self {
        validate_plot_options(Self::default().resolution_points, level);
        Self {
            align_level: level,
            ..Self::default()
        }
    }

    pub(crate) fn validate(&self) {
        validate_plot_options(self.resolution_points, self.align_level);
    }
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
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum DistanceMetric {
    Indel,
    Edit,
    GapLinear,
    GapAffine,
    GapAffine2p,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
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

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum WfaError {
    MissingPenaltyModel,
    InvalidPenalties {
        penalties: Penalties,
        reason: &'static str,
    },
    IncompatibleHeuristics {
        distance_metric: DistanceMetric,
        reason: &'static str,
    },
}

impl fmt::Display for WfaError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            WfaError::MissingPenaltyModel => {
                write!(f, "must set a penalty model before building the aligner")
            }
            WfaError::InvalidPenalties { penalties, reason } => {
                write!(f, "invalid penalties {penalties:?}: {reason}")
            }
            WfaError::IncompatibleHeuristics {
                distance_metric,
                reason,
            } => write!(
                f,
                "heuristics are incompatible with {distance_metric:?}: {reason}"
            ),
        }
    }
}

impl std::error::Error for WfaError {}

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
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
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
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
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

pub(crate) fn validate_max_alignment_steps(max_alignment_steps: i32) {
    assert!(
        max_alignment_steps > 0,
        "max_alignment_steps must be positive"
    );
}

pub(crate) fn validate_max_memory(max_memory_resident: u64, max_memory_abort: u64) {
    assert!(
        max_memory_resident <= max_memory_abort,
        "max_memory_resident must be less than or equal to max_memory_abort"
    );
}

pub(crate) fn validate_max_num_threads(max_num_threads: i32) {
    assert!(max_num_threads > 0, "max_num_threads must be positive");
}

pub(crate) fn validate_min_offsets_per_thread(min_offsets_per_thread: i32) {
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

fn validate_plot_options(resolution_points: i32, align_level: i32) {
    assert!(resolution_points > 0, "resolution_points must be positive");
    assert!(
        align_level >= -1,
        "align_level must be greater than or equal to -1"
    );
}

fn check_i32(value: i64) -> bool {
    value >= i32::MIN as i64 && value <= i32::MAX as i64
}

fn check_adjusted_penalties_fit(penalties: Penalties) -> Result<(), WfaError> {
    let fits = match penalties {
        Penalties::Indel | Penalties::Edit => true,
        Penalties::Linear {
            match_,
            mismatch,
            indel,
        } if match_ < 0 => {
            check_i32(2 * mismatch as i64 - 2 * match_ as i64)
                && check_i32(2 * indel as i64 - match_ as i64)
        }
        Penalties::Affine {
            match_,
            mismatch,
            gap_opening,
            gap_extension,
        } if match_ < 0 => {
            check_i32(2 * mismatch as i64 - 2 * match_ as i64)
                && check_i32(2 * gap_opening as i64)
                && check_i32(2 * gap_extension as i64 - match_ as i64)
        }
        Penalties::Affine2p {
            match_,
            mismatch,
            gap_opening1,
            gap_extension1,
            gap_opening2,
            gap_extension2,
        } if match_ < 0 => {
            check_i32(2 * mismatch as i64 - 2 * match_ as i64)
                && check_i32(2 * gap_opening1 as i64)
                && check_i32(2 * gap_extension1 as i64 - match_ as i64)
                && check_i32(2 * gap_opening2 as i64)
                && check_i32(2 * gap_extension2 as i64 - match_ as i64)
        }
        _ => true,
    };

    if fits {
        Ok(())
    } else {
        Err(WfaError::InvalidPenalties {
            penalties,
            reason: "adjusted penalties must fit in i32",
        })
    }
}

pub(crate) fn validate_penalties(penalties: Penalties) -> Result<(), WfaError> {
    let invalid = |reason| WfaError::InvalidPenalties { penalties, reason };

    match penalties {
        Penalties::Indel | Penalties::Edit => Ok(()),
        Penalties::Linear {
            match_,
            mismatch,
            indel,
        } => {
            if match_ > 0 {
                Err(invalid("match score must be negative or zero"))
            } else if mismatch <= 0 || indel <= 0 {
                Err(invalid(
                    "linear penalties require mismatch > 0 and indel > 0",
                ))
            } else {
                check_adjusted_penalties_fit(penalties)
            }
        }
        Penalties::Affine {
            match_,
            mismatch,
            gap_opening,
            gap_extension,
        } => {
            if match_ > 0 {
                Err(invalid("match score must be negative or zero"))
            } else if mismatch <= 0 || gap_opening < 0 || gap_extension <= 0 {
                Err(invalid(
                    "affine penalties require mismatch > 0, gap_opening >= 0, and gap_extension > 0",
                ))
            } else {
                check_adjusted_penalties_fit(penalties)
            }
        }
        Penalties::Affine2p {
            match_,
            mismatch,
            gap_opening1,
            gap_extension1,
            gap_opening2,
            gap_extension2,
        } => {
            if match_ > 0 {
                return Err(invalid("match score must be negative or zero"));
            }

            if match_ < 0 {
                // WFA2 uses two different mismatch formulas, and we mirror both
                // exactly. The positivity check below uses `2 * mismatch - match`
                // (wavefront_penalties.c: `(2*X - M) > 0`), while the i32 overflow
                // check in `check_adjusted_penalties_fit` uses the stored adjusted
                // value `2 * mismatch - 2 * match`. They are intentionally
                // different, do not make them the same!
                let mismatch_adjusted = 2 * mismatch as i64 - match_ as i64;
                let extension1_adjusted = 2 * gap_extension1 as i64 - match_ as i64;
                let extension2_adjusted = 2 * gap_extension2 as i64 - match_ as i64;
                if mismatch_adjusted <= 0
                    || gap_opening1 < 0
                    || extension1_adjusted <= 0
                    || gap_opening2 < 0
                    || extension2_adjusted <= 0
                {
                    Err(invalid(
                        "affine2p penalties with negative match require (2 * mismatch - match) > 0, gap openings >= 0, and (2 * gap_extension - match) > 0",
                    ))
                } else {
                    check_adjusted_penalties_fit(penalties)
                }
            } else if mismatch <= 0
                || gap_opening1 < 0
                || gap_extension1 <= 0
                || gap_opening2 < 0
                || gap_extension2 <= 0
            {
                Err(invalid(
                    "affine2p penalties require mismatch > 0, gap openings >= 0, and gap extensions > 0",
                ))
            } else {
                Ok(())
            }
        }
    }
}

#[cfg(unix)]
pub(crate) fn validate_heuristics_for_distance_metric(
    heuristics: &Heuristics,
    distance_metric: DistanceMetric,
) {
    check_heuristics_for_distance_metric(heuristics, distance_metric).unwrap_or_else(|err| {
        panic!("{err}");
    });
}

pub(crate) fn check_heuristics_for_distance_metric(
    heuristics: &Heuristics,
    distance_metric: DistanceMetric,
) -> Result<(), WfaError> {
    if heuristics.drop_heuristic().is_some()
        && matches!(
            distance_metric,
            DistanceMetric::Indel | DistanceMetric::Edit
        )
    {
        Err(WfaError::IncompatibleHeuristics {
            distance_metric,
            reason: "drop heuristics are not compatible with edit or indel distance metrics",
        })
    } else {
        Ok(())
    }
}
