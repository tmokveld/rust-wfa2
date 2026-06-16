pub mod aligner;

pub mod wfa2 {
    //! Re-export wfa2-sys bindings
    pub use wfa2_sys::*;
}

#[cfg(test)]
mod tests {
    use crate::aligner::{AlignmentScope, AlignmentStatus, MemoryModel, WFAligner};

    /// Compress WFA2 cigar so that it's easier to read
    fn compress_cigar(cigar: &[u8]) -> String {
        let mut out = String::new();

        let mut runlen = 0;
        let mut current_type = 0;

        for c in cigar {
            if current_type == *c {
                runlen += 1;
            } else {
                if runlen > 0 {
                    out += format!(
                        "{}{}",
                        runlen,
                        std::str::from_utf8(&[current_type]).unwrap()
                    )
                    .as_str();
                }
                runlen = 1;
                current_type = *c;
            }
        }
        if runlen > 0 {
            out += format!(
                "{}{}",
                runlen,
                std::str::from_utf8(&[current_type]).unwrap()
            )
            .as_str();
        }
        out
    }

    fn cigar_string(aligner: &WFAligner) -> String {
        String::from_utf8(aligner.cigar_operations()).unwrap()
    }

    /// Reproduce basic test from library README
    #[test]
    fn test_end_to_end() {
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
            .affine(6, 4, 2)
            .build();

        let pattern = b"TCTTTACTCGCGCGTTGGAGAAATACAATAGT";
        let text = b"TCTATACTGCGCGTTTGGAGAAATAAAATAGT";
        let status = aligner.align_end_to_end(pattern, text);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), -24);
        assert_eq!(
            aligner.cigar_operations(),
            b"MMMXMMMMDMMMMMMMIMMMMMMMMMXMMMMMM"
        );
    }

    /// Test align_ends_free method and affine_with_match builder configuration
    #[test]
    fn test_ends_free() {
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
            .affine_with_match(-1, 3, 2, 1)
            .build();

        let pattern = b"CGCGTTTGGAGAA";
        let text = b"TCTATACTGCGCGTTTGGAGAAATAAAATAGT";
        let pattern_size = pattern.len() as i32;
        let text_size = text.len() as i32;
        let status = aligner.align_ends_free(
            pattern,
            pattern_size,
            pattern_size,
            text,
            text_size,
            text_size,
        );
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), 13);

        // CIGAR output is configured for a reversed notion of pattern/text:
        assert_eq!(compress_cigar(&aligner.cigar_operations()), "9I13M10I");

        let pattern = b"TCTATACTGCGCGTTTGGAGAAATAAAATAGT";
        let text = b"CGCGTTTGGAGAA";
        let pattern_size = pattern.len() as i32;
        let text_size = text.len() as i32;
        let status = aligner.align_ends_free(
            pattern,
            pattern_size,
            pattern_size,
            text,
            text_size,
            text_size,
        );
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), 13);
        assert_eq!(compress_cigar(&aligner.cigar_operations()), "9D13M10D");
    }

    /// Change pattern to test test the left-right shift behavior of this library
    #[test]
    fn test_ends_free_shift() {
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
            .affine_with_match(-1, 3, 2, 1)
            .build();

        let pattern = b"TATATTTTTTTTGGAGAAATAAAATA";
        let text = b"TCTATATTTTTTTTTGGAGAAATAAAATAGT";
        let pattern_size = pattern.len() as i32;
        let text_size = text.len() as i32;
        let status = aligner.align_ends_free(
            pattern,
            pattern_size,
            pattern_size,
            text,
            text_size,
            text_size,
        );
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        //assert_eq!(aligner.score(), 18);

        // CIGAR output is configured for a reversed notion of pattern/text:
        assert_eq!(cigar_string(&aligner), "IIMMMMMMMMMMMMIMMMMMMMMMMMMMMII");
    }

    /// Test double affine mode, and with 0 gap extension to see if a long gap
    /// is created as expected
    #[test]
    fn test_end_to_end_affine2() {
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
            .affine2p_with_match(-1, 3, 3, 3, 10, 0)
            .build();

        let pattern = b"TCTATAATAGT";
        let text = b"TCTATACTGCGCGTTTGGAGAAATAAAATAGT";
        let status = aligner.align_end_to_end(pattern, text);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), 1);
        assert_eq!(compress_cigar(&aligner.cigar_operations()), "6M21I5M");
    }

    /// Test double affine mode, and with 0 gap open
    #[test]
    fn test_end_to_end_affine2_zero_open() {
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
            .affine2p_with_match(-1, 3, 0, 4, 0, 10)
            .build();

        let pattern = b"TCTATAATAGT";
        let text = b"TCTATACTGCGCGTTTGGAGAAATAAAATAGT";
        let status = aligner.align_end_to_end(pattern, text);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), -73);
        assert_eq!(cigar_string(&aligner), "MMMMMMIIIIIIIIIIIIMIIIIMMIIIIIMM");
    }

    /// This test reproduces a bug found in WFA2-lib main branch at 94bcccd.
    ///
    /// A version directly in C is submitted here:
    /// https://github.com/smarco/WFA2-lib/issues/73
    ///
    #[test]
    fn test_ends_free_bug() {
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .linear_with_match(-1, 1, 1)
            .build();

        let pattern = b"A";
        let text = b"ACG";
        let status = aligner.align_ends_free(pattern, 0, 0, text, 0, 2);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), 1);

        // bug version output
        //assert_eq!(aligner.score(), 2);

        // CIGAR output is configured for a reversed notion of pattern/text:
        assert_eq!(cigar_string(&aligner), "MII");
    }

    /// Test simple end to end affine gap example to verify that large gap is created as expected
    ///
    #[test]
    fn test_end_to_end_affine() {
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
            .affine_with_match(-1, 2, 2, 1)
            .build();

        let pattern = b"ATAATA";
        let text = b"ATACATAAAATA";
        let status = aligner.align_end_to_end(pattern, text);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), -2);
        assert_eq!(cigar_string(&aligner), "MMMIIIIIIMMM");
    }

    /// Test case expected to have equal score
    #[test]
    fn test_linear() {
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
            .affine_with_match(-1, 2, 0, 1)
            .build();

        let pattern = b"ATAATA";
        let text = b"ATACATAAAATA";
        let status = aligner.align_end_to_end(pattern, text);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), 0);

        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
            .linear_with_match(-1, 2, 1)
            .build();

        let pattern = b"ATAATA";
        let text = b"ATACATAAAATA";
        let status = aligner.align_end_to_end(pattern, text);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), 0);
    }

    /// Test case expected to have equal score
    #[test]
    fn test_score_only() {
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
            .linear_with_match(-1, 2, 1)
            .build();

        let pattern = b"ATAATA";
        let text = b"ATACATAAAATA";
        let status = aligner.align_end_to_end(pattern, text);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), 0);

        let mut aligner = WFAligner::builder(AlignmentScope::Score, MemoryModel::MemoryLow)
            .linear_with_match(-1, 2, 1)
            .build();

        let pattern = b"ATAATA";
        let text = b"ATACATAAAATA";
        let status = aligner.align_end_to_end(pattern, text);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), 0);
    }
}
