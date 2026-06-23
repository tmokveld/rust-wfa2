/// Derives the aligned span (start and end on both axes) from the active CIGAR operations.
///
/// Both leading and trailing indels are stripped, so the span covers only the aligned core
/// (from the first to the last `M`/`X` column). This keeps the two axes symmetric and is used
/// for ends-free/local alignments. The span is derived from the active CIGAR so it remains
/// consistent across both unidirectional and BiWFA paths.
pub(crate) fn alignment_span_from_ops(raw_operations: &[u8]) -> ((usize, usize), (usize, usize)) {
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
pub(crate) fn extension_alignment_span_from_ops(
    raw_operations: &[u8],
) -> ((usize, usize), (usize, usize)) {
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
