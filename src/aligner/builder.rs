use super::attributes::WFAttributes;
use super::config::{
    check_heuristics_for_distance_metric, validate_memory_model_compatibility, AlignmentScope,
    Heuristics, MemoryModel, PlotOptions, ResourceLimits, WfaError,
};
use super::facade::WFAligner;
use super::raw::WfaRawHandle;

pub struct WFAlignerBuilder {
    attributes: WFAttributes,
    penalty_set: bool,
    heuristics: Option<Heuristics>,
}

impl WFAlignerBuilder {
    /// Create a builder for a WFA2 aligner.
    ///
    /// The builder starts from WFA2's C default attributes with one intentional
    /// semantic change: heuristics are disabled by default. This keeps Rust
    /// alignments exact unless callers opt into heuristic pruning with
    /// [`WFAlignerBuilder::with_heuristics`]. Use
    /// [`Heuristics::wfa2_default`] to recover the C default WF-adaptive
    /// heuristic configuration.
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

    /// Enable WFA2's native `.plot` dump recorder.
    ///
    /// Call [`WFAligner::write_plot`] after an alignment run to write the dump
    /// to a file. This is a debugging/tooling format, not an image renderer.
    pub fn with_plotting(mut self, plot_options: PlotOptions) -> Self {
        self.attributes = self.attributes.plotting(plot_options);
        self
    }

    /// Build the WFAligner with the configured settings.
    pub fn build(self) -> Result<WFAligner, WfaError> {
        if !self.penalty_set {
            return Err(WfaError::MissingPenaltyModel);
        }

        if let Some(heuristics) = self.heuristics {
            check_heuristics_for_distance_metric(&heuristics, self.attributes.distance_metric())?;
        }

        let heuristics = self.heuristics.unwrap_or_default();
        validate_memory_model_compatibility(
            self.attributes.selected_memory_model(),
            self.attributes.alignment_scope_value(),
            self.attributes.penalties(),
            &heuristics,
        )?;

        let mut raw = WfaRawHandle::new(self.attributes)?;

        if let Some(heuristics) = self.heuristics {
            raw.set_heuristics(heuristics);
        }

        Ok(WFAligner { raw })
    }
}
