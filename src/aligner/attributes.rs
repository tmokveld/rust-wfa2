use crate::wfa2;

use super::config::{
    validate_max_alignment_steps, validate_max_memory, validate_max_num_threads,
    validate_min_offsets_per_thread, AlignmentScope, DistanceMetric, MemoryModel, Penalties,
    PlotOptions, ResourceLimits,
};

#[derive(Debug, Copy, Clone)]
pub(crate) struct WFAttributes {
    pub(crate) inner: wfa2::wavefront_aligner_attr_t,
}

impl WFAttributes {
    pub(crate) fn default() -> Self {
        let mut inner = unsafe { wfa2::wavefront_aligner_attr_default };
        // The C default enables WF-adaptive heuristics. The Rust safe wrapper
        // deliberately defaults to exact alignment; callers can opt back into
        // the C default with `Heuristics::wfa2_default()`.
        inner.heuristic.strategy = wfa2::wf_heuristic_strategy_wf_heuristic_none;
        Self { inner }
    }

    pub(crate) fn memory_model(mut self, memory_model: MemoryModel) -> Self {
        let memory_mode = match memory_model {
            MemoryModel::MemoryHigh => wfa2::wavefront_memory_t_wavefront_memory_high,
            MemoryModel::MemoryMed => wfa2::wavefront_memory_t_wavefront_memory_med,
            MemoryModel::MemoryLow => wfa2::wavefront_memory_t_wavefront_memory_low,
            MemoryModel::MemoryUltraLow => wfa2::wavefront_memory_t_wavefront_memory_ultralow,
            MemoryModel::MemorySingletrack => wfa2::wavefront_memory_t_wavefront_memory_singletrack,
        };
        self.inner.memory_mode = memory_mode;
        self
    }

    pub(crate) fn alignment_scope(mut self, alignment_scope: AlignmentScope) -> Self {
        let alignment_scope = match alignment_scope {
            AlignmentScope::Score => wfa2::alignment_scope_t_compute_score,
            AlignmentScope::Alignment => wfa2::alignment_scope_t_compute_alignment,
        };
        self.inner.alignment_scope = alignment_scope;
        self
    }

    pub(crate) fn selected_memory_model(&self) -> MemoryModel {
        match self.inner.memory_mode {
            wfa2::wavefront_memory_t_wavefront_memory_high => MemoryModel::MemoryHigh,
            wfa2::wavefront_memory_t_wavefront_memory_med => MemoryModel::MemoryMed,
            wfa2::wavefront_memory_t_wavefront_memory_low => MemoryModel::MemoryLow,
            wfa2::wavefront_memory_t_wavefront_memory_ultralow => MemoryModel::MemoryUltraLow,
            wfa2::wavefront_memory_t_wavefront_memory_singletrack => MemoryModel::MemorySingletrack,
            _ => panic!("Unknown memory model: {}", self.inner.memory_mode),
        }
    }

    pub(crate) fn alignment_scope_value(&self) -> AlignmentScope {
        AlignmentScope::from(self.inner.alignment_scope)
    }

    pub(crate) fn resource_limits(&self) -> ResourceLimits {
        let system = &self.inner.system;
        ResourceLimits {
            max_alignment_steps: system.max_alignment_steps,
            max_memory_resident: system.max_memory_resident,
            max_memory_abort: system.max_memory_abort,
            max_num_threads: system.max_num_threads,
            min_offsets_per_thread: system.min_offsets_per_thread,
        }
    }

    pub(crate) fn distance_metric(&self) -> DistanceMetric {
        DistanceMetric::from(self.inner.distance_metric)
    }

    pub(crate) fn penalties(&self) -> Penalties {
        match self.distance_metric() {
            DistanceMetric::Indel => Penalties::Indel,
            DistanceMetric::Edit => Penalties::Edit,
            DistanceMetric::GapLinear => Penalties::Linear {
                match_: self.inner.linear_penalties.match_,
                mismatch: self.inner.linear_penalties.mismatch,
                indel: self.inner.linear_penalties.indel,
            },
            DistanceMetric::GapAffine => Penalties::Affine {
                match_: self.inner.affine_penalties.match_,
                mismatch: self.inner.affine_penalties.mismatch,
                gap_opening: self.inner.affine_penalties.gap_opening,
                gap_extension: self.inner.affine_penalties.gap_extension,
            },
            DistanceMetric::GapAffine2p => Penalties::Affine2p {
                match_: self.inner.affine2p_penalties.match_,
                mismatch: self.inner.affine2p_penalties.mismatch,
                gap_opening1: self.inner.affine2p_penalties.gap_opening1,
                gap_extension1: self.inner.affine2p_penalties.gap_extension1,
                gap_opening2: self.inner.affine2p_penalties.gap_opening2,
                gap_extension2: self.inner.affine2p_penalties.gap_extension2,
            },
        }
    }

    pub(crate) fn set_resource_limits(mut self, resource_limits: ResourceLimits) -> Self {
        self = self.max_alignment_steps(resource_limits.max_alignment_steps);
        self = self.max_memory(
            resource_limits.max_memory_resident,
            resource_limits.max_memory_abort,
        );
        self = self.max_num_threads(resource_limits.max_num_threads);
        self.min_offsets_per_thread(resource_limits.min_offsets_per_thread)
    }

    pub(crate) fn max_alignment_steps(mut self, max_alignment_steps: i32) -> Self {
        validate_max_alignment_steps(max_alignment_steps);
        self.inner.system.max_alignment_steps = max_alignment_steps;
        self
    }

    pub(crate) fn max_memory(mut self, max_memory_resident: u64, max_memory_abort: u64) -> Self {
        validate_max_memory(max_memory_resident, max_memory_abort);
        self.inner.system.max_memory_resident = max_memory_resident;
        self.inner.system.max_memory_abort = max_memory_abort;
        self
    }

    pub(crate) fn max_num_threads(mut self, max_num_threads: i32) -> Self {
        validate_max_num_threads(max_num_threads);
        self.inner.system.max_num_threads = max_num_threads;
        self
    }

    pub(crate) fn min_offsets_per_thread(mut self, min_offsets_per_thread: i32) -> Self {
        validate_min_offsets_per_thread(min_offsets_per_thread);
        self.inner.system.min_offsets_per_thread = min_offsets_per_thread;
        self
    }

    pub(crate) fn plotting(mut self, plot_options: PlotOptions) -> Self {
        plot_options.validate();
        self.inner.plot.enabled = true;
        self.inner.plot.resolution_points = plot_options.resolution_points;
        self.inner.plot.align_level = plot_options.align_level;
        self
    }

    pub(crate) fn indel_penalties(mut self) -> Self {
        self.inner.distance_metric = wfa2::distance_metric_t_indel;
        self
    }

    pub(crate) fn edit_penalties(mut self) -> Self {
        self.inner.distance_metric = wfa2::distance_metric_t_edit;
        self
    }

    pub(crate) fn linear_penalties(mut self, match_: i32, mismatch: i32, indel: i32) -> Self {
        self.inner.distance_metric = wfa2::distance_metric_t_gap_linear;
        self.inner.linear_penalties.match_ = match_; // (Penalty representation usually M <= 0)
        self.inner.linear_penalties.mismatch = mismatch; // (Penalty representation usually X > 0)
        self.inner.linear_penalties.indel = indel; // (Penalty representation usually I > 0)
        self
    }

    pub(crate) fn affine_penalties(
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

    pub(crate) fn affine2p_penalties(
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
