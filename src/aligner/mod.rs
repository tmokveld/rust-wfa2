mod attributes;
mod builder;
mod cigar;
mod config;
mod facade;
mod lambda;
mod packed2bits;
mod raw;
mod span;

#[cfg(test)]
mod tests;

pub use builder::WFAlignerBuilder;
pub use cigar::CigarOp;
pub use config::{
    AdaptiveHeuristic, AlignmentResult, AlignmentScope, AlignmentStatus, BandHeuristic,
    DistanceMetric, DropHeuristic, Heuristics, MemoryModel, Penalties, PlotOptions, ResourceLimits,
    WfaAlign, WfaError, WfaOp,
};
pub use facade::WFAligner;
pub use packed2bits::pack_dna_2bits;

#[cfg(test)]
pub(crate) use cigar::CigarView;
#[cfg(test)]
pub(crate) use span::{alignment_span_from_ops, extension_alignment_span_from_ops};
