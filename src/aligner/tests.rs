use super::*;
use std::any::Any;
use std::fs;
use std::io;
use std::os::raw::c_char;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};
use WfaOp::*;

const PATTERN: &[u8] = b"AGCTAGTGTCAATGGCTACTTTTCAGGTCCT";
const TEXT: &[u8] = b"AACTAAGTGTCGGTGGCTACTATATATCAGGTCCT";
static PLOT_TEST_COUNTER: AtomicU64 = AtomicU64::new(0);

fn raw_cigar_string(aligner: &WFAligner) -> String {
    String::from_utf8(aligner.wfa_cigar_bytes()).unwrap()
}

fn panic_message<'a>(payload: &'a (dyn Any + Send + 'static)) -> Option<&'a str> {
    payload
        .downcast_ref::<&str>()
        .copied()
        .or_else(|| payload.downcast_ref::<String>().map(String::as_str))
}

fn temp_plot_path(test_name: &str) -> PathBuf {
    let id = PLOT_TEST_COUNTER.fetch_add(1, Ordering::Relaxed);
    std::env::temp_dir().join(format!(
        "rust_wfa2_{test_name}_{}_{}.plot",
        std::process::id(),
        id
    ))
}

fn read_plot(path: &Path) -> String {
    let plot = fs::read_to_string(path).unwrap();
    let _ = fs::remove_file(path);
    plot
}

fn assert_plot_metadata_and_heatmap(plot: &str) {
    assert!(plot.contains("# PatternLength "));
    assert!(plot.contains("# Pattern "));
    assert!(plot.contains("# TextLength "));
    assert!(plot.contains("# Text "));
    assert!(plot.contains("# Heatmap M\n"));
}

fn assert_alignment_surfaces_match(expected: &mut WFAligner, actual: &mut WFAligner) {
    assert_eq!(actual.score(), expected.score());
    assert_eq!(actual.wfa_cigar_bytes(), expected.wfa_cigar_bytes());
    assert_eq!(actual.sam_cigar_bytes(), expected.sam_cigar_bytes());
    assert_eq!(
        actual.wfa_packed_cigar(true),
        expected.wfa_packed_cigar(true)
    );
    assert_eq!(
        actual.sam_packed_cigar(true),
        expected.sam_packed_cigar(true)
    );
    assert_eq!(actual.count_matches(), expected.count_matches());
    assert_eq!(
        actual.cigar_score_clipped(0),
        expected.cigar_score_clipped(0)
    );
    assert_eq!(actual.cigar_score(), expected.cigar_score());
    assert_eq!(actual.get_alignment_span(), expected.get_alignment_span());
    assert_eq!(actual.get_alignment(), expected.get_alignment());
}

fn run_invalid_penalty_child(case: &str) {
    let result = match case {
        "linear" => WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .linear(0, 1)
            .build(),
        "affine" => WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .affine(0, -1, 0)
            .build(),
        "affine2p" => WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .affine2p(0, 0, 1, 0, 1)
            .build(),
        _ => panic!("unknown invalid penalty case: {case}"),
    };

    assert!(
        matches!(result, Err(WfaError::InvalidPenalties { .. })),
        "expected invalid penalties error for {case}"
    );
}

#[test]
fn test_invalid_penalties_do_not_exit_process() {
    const CHILD_ENV: &str = "RUST_WFA2_INVALID_PENALTY_CHILD";

    if let Ok(case) = std::env::var(CHILD_ENV) {
        run_invalid_penalty_child(&case);
        return;
    }

    for case in ["linear", "affine", "affine2p"] {
        let output = Command::new(std::env::current_exe().unwrap())
            .arg("invalid_penalties_do_not_exit_process")
            .env(CHILD_ENV, case)
            .output()
            .unwrap();

        assert!(
            output.status.success(),
            "invalid {case} penalties exited the process: status={:?}, stderr={}",
            output.status.code(),
            String::from_utf8_lossy(&output.stderr)
        );
    }
}

#[test]
fn test_reap_preserves_aligner_for_reuse() {
    let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .edit()
        .build()
        .unwrap();

    aligner.reap();

    let first_result = aligner.align_end_to_end(PATTERN, TEXT);
    let first_cigar = aligner.wfa_cigar_bytes();
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
        aligner.wfa_cigar_bytes(),
        b"MMMXMMMMDMMMMMMMIMMMMMMMMMXMMMMMM"
    );
}

#[test]
fn test_write_plot_alignment_scope_contains_metadata_heatmap_and_cigar_lists() {
    let path = temp_plot_path("alignment_scope");
    let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .with_plotting(PlotOptions::default())
        .edit()
        .build()
        .unwrap();

    let result = aligner.align_end_to_end(PATTERN, TEXT);
    assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
    aligner.write_plot(&path).unwrap();

    let plot = read_plot(&path);
    assert_plot_metadata_and_heatmap(&plot);
    assert!(plot.contains("# List CIGAR-M "));
    assert!(plot.contains("# List CIGAR-X "));
    assert!(plot.contains("# List CIGAR-I "));
    assert!(plot.contains("# List CIGAR-D "));
}

#[test]
fn test_write_plot_rejects_disabled_plotting() {
    let path = temp_plot_path("disabled");
    let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .edit()
        .build()
        .unwrap();
    let result = aligner.align_end_to_end(PATTERN, TEXT);
    assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);

    let err = aligner.write_plot(&path).unwrap_err();
    assert_eq!(err.kind(), io::ErrorKind::InvalidInput);
}

#[test]
fn test_write_plot_rejects_missing_alignment_run() {
    let path = temp_plot_path("missing_alignment");
    let aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .with_plotting(PlotOptions::default())
        .edit()
        .build()
        .unwrap();

    let err = aligner.write_plot(&path).unwrap_err();
    assert_eq!(err.kind(), io::ErrorKind::InvalidInput);
}

#[test]
fn test_write_plot_score_scope_omits_cigar_lists() {
    let path = temp_plot_path("score_scope");
    let mut aligner = WFAligner::builder(AlignmentScope::Score, MemoryModel::MemoryHigh)
        .with_plotting(PlotOptions::default())
        .edit()
        .build()
        .unwrap();

    let result = aligner.align_end_to_end(PATTERN, TEXT);
    assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
    aligner.write_plot(&path).unwrap();

    let plot = read_plot(&path);
    assert_plot_metadata_and_heatmap(&plot);
    assert!(!plot.contains("# List CIGAR-M "));
}

#[test]
fn test_write_plot_supports_biwfa() {
    let path = temp_plot_path("biwfa");
    let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryUltraLow)
        .with_plotting(PlotOptions::final_alignment())
        .affine(1, 5, 1)
        .build()
        .unwrap();

    let result = aligner.align_end_to_end(PATTERN, TEXT);
    assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
    aligner.write_plot(&path).unwrap();

    let plot = read_plot(&path);
    assert_plot_metadata_and_heatmap(&plot);
}

#[test]
fn test_aligner_indel() {
    let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .indel()
        .build()
        .unwrap();
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
        .build()
        .unwrap();
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
fn test_pack_dna_2bits_layout_and_input_handling() {
    assert_eq!(pack_dna_2bits(b""), Vec::<u8>::new());
    assert_eq!(pack_dna_2bits(b"ACGT"), vec![0b1110_0100]);
    assert_eq!(pack_dna_2bits(b"ACGTA"), vec![0b1110_0100, 0]);
    assert_eq!(pack_dna_2bits(b"acgt"), vec![0b1110_0100]);

    let result = std::panic::catch_unwind(|| pack_dna_2bits(b"ACGN"));
    let payload = result.expect_err("expected invalid DNA base panic");
    assert_eq!(
        panic_message(payload.as_ref()),
        Some("invalid DNA base for 2-bit packing: 0x4E")
    );
}

#[test]
fn test_align_end_to_end_packed2bits_matches_byte_alignment_outputs() {
    let mut byte_aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .affine(1, 5, 1)
        .build()
        .unwrap();
    let byte_result = byte_aligner.align_end_to_end(PATTERN, TEXT);
    assert_eq!(byte_result.status, AlignmentStatus::StatusAlgCompleted);

    let byte_score = byte_aligner.score();
    let byte_wfa_cigar = byte_aligner.wfa_cigar_bytes();
    let byte_sam_cigar = byte_aligner.sam_cigar_bytes();
    let byte_wfa_packed = byte_aligner.wfa_packed_cigar(true);
    let byte_sam_packed = byte_aligner.sam_packed_cigar(true);
    let byte_match_count = byte_aligner.count_matches();
    let byte_clipped_score = byte_aligner.cigar_score_clipped(0);
    let byte_cigar_score = byte_aligner.cigar_score();
    let byte_span = byte_aligner.get_alignment_span();
    let byte_alignment = byte_aligner.get_alignment();

    let mut packed_pattern = pack_dna_2bits(PATTERN);
    let mut packed_text = pack_dna_2bits(TEXT);
    packed_pattern.push(0xFF);
    packed_text.push(0xFF);

    let mut packed_aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .affine(1, 5, 1)
        .build()
        .unwrap();
    let packed_result = packed_aligner.align_end_to_end_packed2bits(
        &packed_pattern,
        PATTERN.len(),
        &packed_text,
        TEXT.len(),
    );

    assert_eq!(packed_result.status, byte_result.status);
    assert_eq!(packed_result.score, byte_result.score);
    assert_eq!(packed_aligner.score(), byte_score);
    assert_eq!(packed_aligner.wfa_cigar_bytes(), byte_wfa_cigar);
    assert_eq!(packed_aligner.sam_cigar_bytes(), byte_sam_cigar);
    assert_eq!(packed_aligner.wfa_packed_cigar(true), byte_wfa_packed);
    assert_eq!(packed_aligner.sam_packed_cigar(true), byte_sam_packed);
    assert_eq!(packed_aligner.count_matches(), byte_match_count);
    assert_eq!(packed_aligner.cigar_score_clipped(0), byte_clipped_score);
    assert_eq!(packed_aligner.cigar_score(), byte_cigar_score);
    assert_eq!(packed_aligner.get_alignment_span(), byte_span);
    assert_eq!(packed_aligner.get_alignment(), byte_alignment);
}

#[test]
fn test_align_end_to_end_packed2bits_supports_score_scope_ultralow() {
    let mut byte_aligner = WFAligner::builder(AlignmentScope::Score, MemoryModel::MemoryUltraLow)
        .affine(1, 5, 1)
        .build()
        .unwrap();
    let byte_result = byte_aligner.align_end_to_end(PATTERN, TEXT);
    assert_eq!(byte_result.status, AlignmentStatus::StatusAlgCompleted);

    let packed_pattern = pack_dna_2bits(PATTERN);
    let packed_text = pack_dna_2bits(TEXT);
    let mut packed_aligner = WFAligner::builder(AlignmentScope::Score, MemoryModel::MemoryUltraLow)
        .affine(1, 5, 1)
        .build()
        .unwrap();
    let packed_result = packed_aligner.align_end_to_end_packed2bits(
        &packed_pattern,
        PATTERN.len(),
        &packed_text,
        TEXT.len(),
    );

    assert_eq!(packed_result.status, byte_result.status);
    assert_eq!(packed_result.score, byte_result.score);
    assert_eq!(packed_aligner.score(), byte_aligner.score());
    assert_eq!(packed_aligner.cigar_string(None), "");
}

#[test]
fn test_align_end_to_end_packed2bits_handles_empty_sequences() {
    let text = b"ACGT";
    let packed_empty = pack_dna_2bits(b"");
    let packed_text = pack_dna_2bits(text);

    let mut empty_pattern = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .affine(6, 4, 2)
        .build()
        .unwrap();
    let result =
        empty_pattern.align_end_to_end_packed2bits(&packed_empty, 0, &packed_text, text.len());
    assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
    assert_eq!(raw_cigar_string(&empty_pattern), "IIII");
    let alignment = empty_pattern.get_alignment();
    assert_eq!(alignment.xlen, 0);
    assert_eq!(alignment.ylen, text.len());
    assert_eq!(empty_pattern.get_alignment_span(), ((0, 0), (0, 4)));

    let mut both_empty = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .affine(6, 4, 2)
        .build()
        .unwrap();
    let result = both_empty.align_end_to_end_packed2bits(&packed_empty, 0, &packed_empty, 0);
    assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
    assert_eq!(both_empty.score(), 0);
    assert_eq!(both_empty.get_alignment_span(), ((0, 0), (0, 0)));
    assert!(both_empty.wfa_cigar_bytes().is_empty());
}

#[test]
fn test_align_end_to_end_packed2bits_ignores_unused_tail_bits() {
    let pattern = b"ACGTA";
    let text = b"ACGTA";

    let mut byte_aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .edit()
        .build()
        .unwrap();
    let byte_result = byte_aligner.align_end_to_end(pattern, text);
    assert_eq!(byte_result.status, AlignmentStatus::StatusAlgCompleted);

    let mut packed_pattern = pack_dna_2bits(pattern);
    let mut packed_text = pack_dna_2bits(text);
    packed_pattern[1] = 0b1111_1100;
    packed_text[1] = 0b1111_1100;

    let mut packed_aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .edit()
        .build()
        .unwrap();
    let packed_result = packed_aligner.align_end_to_end_packed2bits(
        &packed_pattern,
        pattern.len(),
        &packed_text,
        text.len(),
    );

    assert_eq!(packed_result.status, byte_result.status);
    assert_eq!(packed_result.score, byte_result.score);
    assert_eq!(
        packed_aligner.wfa_cigar_bytes(),
        byte_aligner.wfa_cigar_bytes()
    );
    assert_eq!(
        packed_aligner.get_alignment_span(),
        byte_aligner.get_alignment_span()
    );
}

#[test]
fn test_align_end_to_end_packed2bits_rejects_too_large_logical_length() {
    let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .edit()
        .build()
        .unwrap();

    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        aligner.align_end_to_end_packed2bits(&[], i32::MAX as usize + 1, &[], 0);
    }));

    let payload = result.expect_err("expected packed logical length panic");
    assert_eq!(
        panic_message(payload.as_ref()),
        Some("pattern logical length must fit in i32")
    );
}

#[test]
fn test_align_end_to_end_packed2bits_rejects_short_backing_slice() {
    let packed_text = pack_dna_2bits(b"A");
    let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .edit()
        .build()
        .unwrap();

    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        aligner.align_end_to_end_packed2bits(&[], 1, &packed_text, 1);
    }));

    let payload = result.expect_err("expected short packed buffer panic");
    assert_eq!(
        panic_message(payload.as_ref()),
        Some("pattern packed 2-bit buffer is too short for logical length")
    );
}

#[test]
fn test_align_end_to_end_lambda_matches_byte_alignment_outputs() {
    let mut byte_aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .affine(1, 5, 1)
        .build()
        .unwrap();
    let byte_result = byte_aligner.align_end_to_end(PATTERN, TEXT);
    assert_eq!(byte_result.status, AlignmentStatus::StatusAlgCompleted);

    let byte_score = byte_aligner.score();
    let byte_wfa_cigar = byte_aligner.wfa_cigar_bytes();
    let byte_sam_cigar = byte_aligner.sam_cigar_bytes();
    let byte_wfa_packed = byte_aligner.wfa_packed_cigar(true);
    let byte_sam_packed = byte_aligner.sam_packed_cigar(true);
    let byte_match_count = byte_aligner.count_matches();
    let byte_clipped_score = byte_aligner.cigar_score_clipped(0);
    let byte_cigar_score = byte_aligner.cigar_score();
    let byte_span = byte_aligner.get_alignment_span();
    let byte_alignment = byte_aligner.get_alignment();

    let matcher_calls = AtomicU64::new(0);
    let mut lambda_aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .affine(1, 5, 1)
        .build()
        .unwrap();
    let lambda_result = lambda_aligner.align_end_to_end_lambda(
        PATTERN.len(),
        TEXT.len(),
        |pattern_pos, text_pos| {
            matcher_calls.fetch_add(1, Ordering::Relaxed);
            PATTERN[pattern_pos] == TEXT[text_pos]
        },
    );

    assert!(matcher_calls.load(Ordering::Relaxed) > 0);
    assert_eq!(lambda_result.status, byte_result.status);
    assert_eq!(lambda_result.score, byte_result.score);
    assert_eq!(lambda_aligner.score(), byte_score);
    assert_eq!(lambda_aligner.wfa_cigar_bytes(), byte_wfa_cigar);
    assert_eq!(lambda_aligner.sam_cigar_bytes(), byte_sam_cigar);
    assert_eq!(lambda_aligner.wfa_packed_cigar(true), byte_wfa_packed);
    assert_eq!(lambda_aligner.sam_packed_cigar(true), byte_sam_packed);
    assert_eq!(lambda_aligner.count_matches(), byte_match_count);
    assert_eq!(lambda_aligner.cigar_score_clipped(0), byte_clipped_score);
    assert_eq!(lambda_aligner.cigar_score(), byte_cigar_score);
    assert_eq!(lambda_aligner.get_alignment_span(), byte_span);
    assert_eq!(lambda_aligner.get_alignment(), byte_alignment);
}

#[test]
fn test_align_end_to_end_lambda_supports_score_scope() {
    let mut byte_aligner = WFAligner::builder(AlignmentScope::Score, MemoryModel::MemoryLow)
        .affine(6, 4, 2)
        .build()
        .unwrap();
    let byte_result = byte_aligner.align_end_to_end(PATTERN, TEXT);
    assert_eq!(byte_result.status, AlignmentStatus::StatusAlgCompleted);

    let mut lambda_aligner = WFAligner::builder(AlignmentScope::Score, MemoryModel::MemoryLow)
        .affine(6, 4, 2)
        .build()
        .unwrap();
    let lambda_result = lambda_aligner.align_end_to_end_lambda(
        PATTERN.len(),
        TEXT.len(),
        |pattern_pos, text_pos| PATTERN[pattern_pos] == TEXT[text_pos],
    );

    assert_eq!(lambda_result.status, byte_result.status);
    assert_eq!(lambda_result.score, byte_result.score);
    assert_eq!(lambda_aligner.score(), byte_aligner.score());
    assert_eq!(lambda_aligner.cigar_string(None), "");
}

#[test]
fn test_align_end_to_end_lambda_resumes_matcher_panic() {
    let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .edit()
        .build()
        .unwrap();

    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        aligner.align_end_to_end_lambda(4, 4, |_, _| -> bool {
            panic!("lambda matcher panic");
        });
    }));

    let payload = result.expect_err("expected matcher panic to resume after WFA2 returns");
    let message = payload
        .downcast_ref::<&str>()
        .copied()
        .or_else(|| payload.downcast_ref::<String>().map(String::as_str));
    assert_eq!(message, Some("lambda matcher panic"));
}

#[test]
fn test_align_end_to_end_lambda_ultralow_reports_original_lengths() {
    let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryUltraLow)
        .affine(1, 5, 1)
        .build()
        .unwrap();

    let result =
        aligner.align_end_to_end_lambda(PATTERN.len(), TEXT.len(), |pattern_pos, text_pos| {
            PATTERN[pattern_pos] == TEXT[text_pos]
        });
    assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);

    let alignment = aligner.get_alignment();
    assert_eq!(alignment.xlen, PATTERN.len());
    assert_eq!(alignment.ylen, TEXT.len());
    assert_eq!(alignment.xstart, 0);
    assert_eq!(alignment.xend, PATTERN.len());
    assert_eq!(alignment.ystart, 0);
    assert_eq!(alignment.yend, TEXT.len());
    assert_eq!(
        aligner.get_alignment_span(),
        ((0, PATTERN.len()), (0, TEXT.len()))
    );
}

#[test]
fn test_aligner_gap_linear() {
    let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .linear(6, 2)
        .build()
        .unwrap();
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
        .build()
        .unwrap();
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
        .build()
        .unwrap();

    let pattern = b"TCTTTACTCGCGCGTTGGAGAAATACAATAGT";
    let text = b"TCTATACTGCGCGTTTGGAGAAATAAAATAGT";
    let result = aligner.align_end_to_end(pattern, text);
    assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
    assert_eq!(aligner.score(), -24);
    assert_eq!(
        aligner.wfa_cigar_bytes(),
        b"MMMXMMMMDMMMMMMMIMMMMMMMMMXMMMMMM"
    );
}

#[test]
fn test_affine_with_match_long_gap() {
    let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
        .affine_with_match(-1, 2, 2, 1)
        .build()
        .unwrap();

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
        .build()
        .unwrap();
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
        .build()
        .unwrap();
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
        .build()
        .unwrap();

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
        .build()
        .unwrap();

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

    let mut affine_aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
        .affine_with_match(-1, 2, 0, 1)
        .build()
        .unwrap();
    let result = affine_aligner.align_end_to_end(pattern, text);
    assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
    assert_eq!(affine_aligner.score(), 0);

    let mut linear_aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
        .linear_with_match(-1, 2, 1)
        .build()
        .unwrap();
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
        .build()
        .unwrap();
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
        .build()
        .unwrap();
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
        .build()
        .unwrap();
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
fn test_align_ends_free_lambda_matches_byte_alignment_outputs() {
    let pattern = b"AATTTAAGTCTAGGCTACTTTC";
    let text = b"CCGACTACTACGAAATTTAAGTATAGGCTACTTTCCGTACGTACGTACGT";
    let text_end_free = text.len() as i32;

    let mut byte_aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .affine(6, 4, 2)
        .build()
        .unwrap();
    let byte_result = byte_aligner.align_ends_free(pattern, 0, 0, text, 0, text_end_free);
    assert_eq!(byte_result.status, AlignmentStatus::StatusAlgCompleted);

    let byte_score = byte_aligner.score();
    let byte_wfa_cigar = byte_aligner.wfa_cigar_bytes();
    let byte_sam_cigar = byte_aligner.sam_cigar_bytes();
    let byte_wfa_packed = byte_aligner.wfa_packed_cigar(true);
    let byte_sam_packed = byte_aligner.sam_packed_cigar(true);
    let byte_span = byte_aligner.get_alignment_span();
    let byte_alignment = byte_aligner.get_alignment();

    let mut lambda_aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .affine(6, 4, 2)
        .build()
        .unwrap();
    let lambda_result = lambda_aligner.align_ends_free_lambda(
        pattern.len(),
        0,
        0,
        text.len(),
        0,
        text_end_free,
        |pattern_pos, text_pos| pattern[pattern_pos] == text[text_pos],
    );

    assert_eq!(lambda_result.status, byte_result.status);
    assert_eq!(lambda_result.score, byte_result.score);
    assert_eq!(lambda_aligner.score(), byte_score);
    assert_eq!(lambda_aligner.wfa_cigar_bytes(), byte_wfa_cigar);
    assert_eq!(lambda_aligner.sam_cigar_bytes(), byte_sam_cigar);
    assert_eq!(lambda_aligner.wfa_packed_cigar(true), byte_wfa_packed);
    assert_eq!(lambda_aligner.sam_packed_cigar(true), byte_sam_packed);
    assert_eq!(lambda_aligner.get_alignment_span(), byte_span);
    assert_eq!(lambda_aligner.get_alignment(), byte_alignment);
}

#[test]
fn test_align_ends_free_lambda_rejects_ultralow_nonzero_free_ends() {
    let matcher_calls = AtomicU64::new(0);
    let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryUltraLow)
        .edit()
        .build()
        .unwrap();

    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        aligner.align_ends_free_lambda(4, 0, 1, 4, 0, 0, |_, _| {
            matcher_calls.fetch_add(1, Ordering::Relaxed);
            true
        });
    }));

    let payload = result.expect_err("expected MemoryUltraLow ends-free lambda rejection");
    let message = payload
        .downcast_ref::<&str>()
        .copied()
        .or_else(|| payload.downcast_ref::<String>().map(String::as_str));
    assert_eq!(
        message,
        Some("Ends-free alignment is not supported with MemoryUltraLow")
    );
    assert_eq!(matcher_calls.load(Ordering::Relaxed), 0);
}

#[test]
fn test_align_ends_free_lambda_ultralow_all_zero_matches_global_lambda() {
    let mut global = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryUltraLow)
        .affine(1, 5, 1)
        .build()
        .unwrap();
    let global_result =
        global.align_end_to_end_lambda(PATTERN.len(), TEXT.len(), |pattern_pos, text_pos| {
            PATTERN[pattern_pos] == TEXT[text_pos]
        });
    assert_eq!(global_result.status, AlignmentStatus::StatusAlgCompleted);
    let global_score = global.score();

    let mut ends_free = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryUltraLow)
        .affine(1, 5, 1)
        .build()
        .unwrap();
    let ends_free_result = ends_free.align_ends_free_lambda(
        PATTERN.len(),
        0,
        0,
        TEXT.len(),
        0,
        0,
        |pattern_pos, text_pos| PATTERN[pattern_pos] == TEXT[text_pos],
    );
    assert_eq!(ends_free_result.status, AlignmentStatus::StatusAlgCompleted);
    assert_eq!(ends_free.score(), global_score);
}

#[test]
fn test_align_ends_free_packed2bits_matches_byte_alignment_outputs() {
    let pattern = b"AATTTAAGTCTAGGCTACTTTC";
    let text = b"CCGACTACTACGAAATTTAAGTATAGGCTACTTTCCGTACGTACGTACGT";
    let text_end_free = text.len() as i32;

    let mut byte_aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .affine(6, 4, 2)
        .build()
        .unwrap();
    let byte_result = byte_aligner.align_ends_free(pattern, 0, 0, text, 0, text_end_free);
    assert_eq!(byte_result.status, AlignmentStatus::StatusAlgCompleted);

    let byte_score = byte_aligner.score();
    let byte_wfa_cigar = byte_aligner.wfa_cigar_bytes();
    let byte_sam_cigar = byte_aligner.sam_cigar_bytes();
    let byte_wfa_packed = byte_aligner.wfa_packed_cigar(true);
    let byte_sam_packed = byte_aligner.sam_packed_cigar(true);
    let byte_span = byte_aligner.get_alignment_span();
    let byte_alignment = byte_aligner.get_alignment();

    let packed_pattern = pack_dna_2bits(pattern);
    let packed_text = pack_dna_2bits(text);
    let mut packed_aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .affine(6, 4, 2)
        .build()
        .unwrap();
    let packed_result = packed_aligner.align_ends_free_packed2bits(
        &packed_pattern,
        pattern.len(),
        0,
        0,
        &packed_text,
        text.len(),
        0,
        text_end_free,
    );

    assert_eq!(packed_result.status, byte_result.status);
    assert_eq!(packed_result.score, byte_result.score);
    assert_eq!(packed_aligner.score(), byte_score);
    assert_eq!(packed_aligner.wfa_cigar_bytes(), byte_wfa_cigar);
    assert_eq!(packed_aligner.sam_cigar_bytes(), byte_sam_cigar);
    assert_eq!(packed_aligner.wfa_packed_cigar(true), byte_wfa_packed);
    assert_eq!(packed_aligner.sam_packed_cigar(true), byte_sam_packed);
    assert_eq!(packed_aligner.get_alignment_span(), byte_span);
    assert_eq!(packed_aligner.get_alignment(), byte_alignment);
}

#[test]
fn test_align_ends_free_packed2bits_rejects_ultralow_nonzero_free_ends() {
    let packed = pack_dna_2bits(b"ACGT");
    let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryUltraLow)
        .edit()
        .build()
        .unwrap();

    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        aligner.align_ends_free_packed2bits(&packed, 4, 0, 1, &packed, 4, 0, 0);
    }));

    let payload = result.expect_err("expected MemoryUltraLow ends-free packed rejection");
    assert_eq!(
        panic_message(payload.as_ref()),
        Some("Ends-free alignment is not supported with MemoryUltraLow")
    );
}

#[test]
fn test_align_ends_free_packed2bits_ultralow_all_zero_matches_global() {
    let packed_pattern = pack_dna_2bits(PATTERN);
    let packed_text = pack_dna_2bits(TEXT);

    let mut global = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryUltraLow)
        .affine(1, 5, 1)
        .build()
        .unwrap();
    let global_result = global.align_end_to_end_packed2bits(
        &packed_pattern,
        PATTERN.len(),
        &packed_text,
        TEXT.len(),
    );
    assert_eq!(global_result.status, AlignmentStatus::StatusAlgCompleted);
    let global_score = global.score();

    let mut ends_free = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryUltraLow)
        .affine(1, 5, 1)
        .build()
        .unwrap();
    let ends_free_result = ends_free.align_ends_free_packed2bits(
        &packed_pattern,
        PATTERN.len(),
        0,
        0,
        &packed_text,
        TEXT.len(),
        0,
        0,
    );
    assert_eq!(ends_free_result.status, AlignmentStatus::StatusAlgCompleted);
    assert_eq!(ends_free.score(), global_score);
}

#[test]
fn test_ends_free_with_match_penalties() {
    let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
        .affine_with_match(-1, 3, 2, 1)
        .build()
        .unwrap();

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
        .build()
        .unwrap();

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
        .build()
        .unwrap();
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
            .build()
            .unwrap();
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
            .build()
            .unwrap();
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
        .build()
        .unwrap();

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
fn test_align_extension_lambda_matches_byte_alignment_outputs() {
    let pattern = b"AATTTAAGTCTGCTACTTTCACGCAGCT";
    let text = b"AATTTCAGTCTGGCTACTTTCACGTACGATGACAGACTCT";

    let mut byte_aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .affine(6, 4, 2)
        .build()
        .unwrap();
    let byte_result = byte_aligner.align_extension(pattern, text);
    assert_eq!(byte_result.status, AlignmentStatus::StatusAlgPartial);

    let byte_score = byte_aligner.score();
    let byte_wfa_cigar = byte_aligner.wfa_cigar_bytes();
    let byte_sam_cigar = byte_aligner.sam_cigar_bytes();
    let byte_wfa_packed = byte_aligner.wfa_packed_cigar(true);
    let byte_sam_packed = byte_aligner.sam_packed_cigar(true);
    let byte_cigar_score = byte_aligner.cigar_score();
    let byte_span = byte_aligner.get_alignment_span();
    let byte_alignment = byte_aligner.get_alignment();

    let mut lambda_aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .affine(6, 4, 2)
        .build()
        .unwrap();
    let lambda_result = lambda_aligner.align_extension_lambda(
        pattern.len(),
        text.len(),
        |pattern_pos, text_pos| pattern[pattern_pos] == text[text_pos],
    );

    assert_eq!(lambda_result.status, byte_result.status);
    assert_eq!(lambda_result.score, byte_result.score);
    assert_eq!(lambda_aligner.score(), byte_score);
    assert_eq!(lambda_aligner.wfa_cigar_bytes(), byte_wfa_cigar);
    assert_eq!(lambda_aligner.sam_cigar_bytes(), byte_sam_cigar);
    assert_eq!(lambda_aligner.wfa_packed_cigar(true), byte_wfa_packed);
    assert_eq!(lambda_aligner.sam_packed_cigar(true), byte_sam_packed);
    assert_eq!(lambda_aligner.cigar_score(), byte_cigar_score);
    assert_eq!(lambda_aligner.get_alignment_span(), byte_span);
    assert_eq!(lambda_aligner.get_alignment(), byte_alignment);
}

#[test]
fn test_align_extension_packed2bits_matches_byte_alignment_outputs() {
    let pattern = b"AATTTAAGTCTGCTACTTTCACGCAGCT";
    let text = b"AATTTCAGTCTGGCTACTTTCACGTACGATGACAGACTCT";

    let mut byte_aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .affine(6, 4, 2)
        .build()
        .unwrap();
    let byte_result = byte_aligner.align_extension(pattern, text);
    assert_eq!(byte_result.status, AlignmentStatus::StatusAlgPartial);

    let byte_score = byte_aligner.score();
    let byte_wfa_cigar = byte_aligner.wfa_cigar_bytes();
    let byte_sam_cigar = byte_aligner.sam_cigar_bytes();
    let byte_wfa_packed = byte_aligner.wfa_packed_cigar(true);
    let byte_sam_packed = byte_aligner.sam_packed_cigar(true);
    let byte_cigar_score = byte_aligner.cigar_score();
    let byte_span = byte_aligner.get_alignment_span();
    let byte_alignment = byte_aligner.get_alignment();

    let packed_pattern = pack_dna_2bits(pattern);
    let packed_text = pack_dna_2bits(text);
    let mut packed_aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .affine(6, 4, 2)
        .build()
        .unwrap();
    let packed_result = packed_aligner.align_extension_packed2bits(
        &packed_pattern,
        pattern.len(),
        &packed_text,
        text.len(),
    );

    assert_eq!(packed_result.status, byte_result.status);
    assert_eq!(packed_result.score, byte_result.score);
    assert_eq!(packed_aligner.score(), byte_score);
    assert_eq!(packed_aligner.wfa_cigar_bytes(), byte_wfa_cigar);
    assert_eq!(packed_aligner.sam_cigar_bytes(), byte_sam_cigar);
    assert_eq!(packed_aligner.wfa_packed_cigar(true), byte_wfa_packed);
    assert_eq!(packed_aligner.sam_packed_cigar(true), byte_sam_packed);
    assert_eq!(packed_aligner.cigar_score(), byte_cigar_score);
    assert_eq!(packed_aligner.get_alignment_span(), byte_span);
    assert_eq!(packed_aligner.get_alignment(), byte_alignment);
}

#[test]
fn test_align_extension_lambda_empty_prefix_has_zero_span() {
    let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .affine(6, 4, 2)
        .build()
        .unwrap();

    let result = aligner.align_extension_lambda(8, 8, |_, _| false);
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
fn test_align_extension_packed2bits_empty_prefix_has_zero_span() {
    let packed_pattern = pack_dna_2bits(b"AAAAAAAA");
    let packed_text = pack_dna_2bits(b"TTTTTTTT");
    let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .affine(6, 4, 2)
        .build()
        .unwrap();

    let result = aligner.align_extension_packed2bits(&packed_pattern, 8, &packed_text, 8);
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
fn test_align_extension_lambda_rejects_ultralow_memory() {
    let matcher_calls = AtomicU64::new(0);
    let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryUltraLow)
        .affine(6, 4, 2)
        .build()
        .unwrap();

    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        aligner.align_extension_lambda(4, 4, |_, _| {
            matcher_calls.fetch_add(1, Ordering::Relaxed);
            true
        });
    }));

    let payload = result.expect_err("expected MemoryUltraLow extension lambda rejection");
    let message = payload
        .downcast_ref::<&str>()
        .copied()
        .or_else(|| payload.downcast_ref::<String>().map(String::as_str));
    assert_eq!(
        message,
        Some("Extension alignment is not supported with MemoryUltraLow")
    );
    assert_eq!(matcher_calls.load(Ordering::Relaxed), 0);
}

#[test]
fn test_align_extension_packed2bits_rejects_ultralow_memory() {
    let packed = pack_dna_2bits(b"ACGT");
    let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryUltraLow)
        .affine(6, 4, 2)
        .build()
        .unwrap();

    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        aligner.align_extension_packed2bits(&packed, 4, &packed, 4);
    }));

    let payload = result.expect_err("expected MemoryUltraLow extension packed rejection");
    assert_eq!(
        panic_message(payload.as_ref()),
        Some("Extension alignment is not supported with MemoryUltraLow")
    );
}

#[test]
fn test_extension_alignment_span_from_ops() {
    // Anchored at the origin: leading and trailing indels extend the span (unlike the
    // ends-free/local span, which strips them).
    assert_eq!(extension_alignment_span_from_ops(b""), ((0, 0), (0, 0)));
    assert_eq!(extension_alignment_span_from_ops(b"MMM"), ((0, 3), (0, 3)));
    assert_eq!(extension_alignment_span_from_ops(b"IMMM"), ((0, 3), (0, 4)));
    assert_eq!(
        extension_alignment_span_from_ops(b"MMMII"),
        ((0, 3), (0, 5))
    );
    assert_eq!(
        extension_alignment_span_from_ops(b"DDMMX"),
        ((0, 5), (0, 3))
    );
}

#[test]
fn test_aligner_extension_supports_score_scope() {
    let pattern = b"AATTTAAGTCTGCTACTTTCACGCAGCT";
    let text = b"AATTTCAGTCTGGCTACTTTCACGTACGATGACAGACTCT";
    let mut aligner = WFAligner::builder(AlignmentScope::Score, MemoryModel::MemoryHigh)
        .affine(6, 4, 2)
        .build()
        .unwrap();

    let result = aligner.align_extension(pattern, text);
    assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
    assert_eq!(aligner.score(), -24);
}

#[test]
#[should_panic(expected = "Extension alignment is not supported with MemoryUltraLow")]
fn test_aligner_extension_rejects_ultralow_memory() {
    let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryUltraLow)
        .affine(6, 4, 2)
        .build()
        .unwrap();

    aligner.align_extension(b"ACGT", b"ACGT");
}

#[test]
#[should_panic(expected = "Ends-free alignment is not supported with MemoryUltraLow")]
fn test_aligner_ends_free_rejects_ultralow_memory() {
    let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryUltraLow)
        .affine(6, 4, 2)
        .build()
        .unwrap();

    aligner.align_ends_free(b"ACGT", 0, 0, b"ACGT", 0, 1);
}

#[test]
fn test_aligner_ends_free_left_extent() {
    let pattern = b"CTTTCACGTACGTGACAGTCTCT";
    let text = b"AATTTCAGTCTGGCTACTTTCACGTACGATGACAGACTCT";
    let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .affine(6, 4, 2)
        .build()
        .unwrap();
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
        .build()
        .unwrap();
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
        .build()
        .unwrap();

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
        .build()
        .unwrap();
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
        .build()
        .unwrap();
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
            .build()
            .unwrap();
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
fn test_singletrack_gap_affine_end_to_end_matches_high_memory() {
    let mut high = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .affine(6, 4, 2)
        .build()
        .unwrap();
    let high_result = high.align_end_to_end(PATTERN, TEXT);
    assert_eq!(high_result.status, AlignmentStatus::StatusAlgCompleted);

    let mut singletrack =
        WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemorySingletrack)
            .affine(6, 4, 2)
            .build()
            .unwrap();
    let singletrack_result = singletrack.align_end_to_end(PATTERN, TEXT);

    assert_eq!(singletrack_result.status, high_result.status);
    assert_eq!(singletrack_result.score, high_result.score);
    assert_alignment_surfaces_match(&mut high, &mut singletrack);
}

#[test]
fn test_singletrack_gap_affine2p_end_to_end_matches_high_memory() {
    let mut high = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .affine2p_with_match(-1, 3, 3, 3, 10, 0)
        .build()
        .unwrap();
    let high_result = high.align_end_to_end(PATTERN, TEXT);
    assert_eq!(high_result.status, AlignmentStatus::StatusAlgCompleted);

    let mut singletrack =
        WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemorySingletrack)
            .affine2p_with_match(-1, 3, 3, 3, 10, 0)
            .build()
            .unwrap();
    let singletrack_result = singletrack.align_end_to_end(PATTERN, TEXT);

    assert_eq!(singletrack_result.status, high_result.status);
    assert_eq!(singletrack_result.score, high_result.score);
    assert_alignment_surfaces_match(&mut high, &mut singletrack);
}

#[test]
fn test_singletrack_ends_free_matches_high_memory() {
    let pattern = b"AGTGTCAATGGCTAC";
    let text = b"GGGGGGGGGGAGTGTCAATGGCTACGGGGGGGGGG";

    let mut high = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .affine_with_match(-1, 2, 2, 1)
        .build()
        .unwrap();
    let high_result =
        high.align_ends_free(pattern, 0, 0, text, text.len() as i32, text.len() as i32);
    assert_eq!(high_result.status, AlignmentStatus::StatusAlgCompleted);

    let mut singletrack =
        WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemorySingletrack)
            .affine_with_match(-1, 2, 2, 1)
            .build()
            .unwrap();
    let singletrack_result =
        singletrack.align_ends_free(pattern, 0, 0, text, text.len() as i32, text.len() as i32);

    assert_eq!(singletrack_result.status, high_result.status);
    assert_eq!(singletrack_result.score, high_result.score);
    assert_alignment_surfaces_match(&mut high, &mut singletrack);
}

#[test]
fn test_singletrack_extension_matches_high_memory() {
    let pattern = b"AATTTAAGTCTGCTACTTTCACGCAGCT";
    let text = b"AATTTCAGTCTGGCTACTTTCACGTACGATGACAGACTCT";

    let mut high = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .affine(6, 4, 2)
        .build()
        .unwrap();
    let high_result = high.align_extension(pattern, text);
    assert_eq!(high_result.status, AlignmentStatus::StatusAlgPartial);

    let mut singletrack =
        WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemorySingletrack)
            .affine(6, 4, 2)
            .build()
            .unwrap();
    let singletrack_result = singletrack.align_extension(pattern, text);

    assert_eq!(singletrack_result.status, high_result.status);
    assert_eq!(singletrack_result.score, high_result.score);
    assert_alignment_surfaces_match(&mut high, &mut singletrack);
}

#[test]
fn test_singletrack_packed2bits_matches_high_memory() {
    let packed_pattern = pack_dna_2bits(PATTERN);
    let packed_text = pack_dna_2bits(TEXT);

    let mut high = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .affine(1, 5, 1)
        .build()
        .unwrap();
    let high_result =
        high.align_end_to_end_packed2bits(&packed_pattern, PATTERN.len(), &packed_text, TEXT.len());
    assert_eq!(high_result.status, AlignmentStatus::StatusAlgCompleted);

    let mut singletrack =
        WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemorySingletrack)
            .affine(1, 5, 1)
            .build()
            .unwrap();
    let singletrack_result = singletrack.align_end_to_end_packed2bits(
        &packed_pattern,
        PATTERN.len(),
        &packed_text,
        TEXT.len(),
    );

    assert_eq!(singletrack_result.status, high_result.status);
    assert_eq!(singletrack_result.score, high_result.score);
    assert_alignment_surfaces_match(&mut high, &mut singletrack);
}

#[test]
fn test_singletrack_allows_non_banded_heuristics() {
    for heuristics in [
        Heuristics::wf_adaptive(1, 10, 50),
        Heuristics::wf_mash(1, 10, 50),
        Heuristics::xdrop(1, 10),
        Heuristics::zdrop(1, 10),
    ] {
        let mut aligner =
            WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemorySingletrack)
                .affine(6, 4, 2)
                .with_heuristics(heuristics)
                .build()
                .unwrap();
        let result = aligner.align_end_to_end(PATTERN, TEXT);
        assert!(matches!(
            result.status,
            AlignmentStatus::StatusAlgCompleted | AlignmentStatus::StatusAlgPartial
        ));
    }
}

#[test]
fn test_singletrack_rejects_unsupported_builder_configs() {
    let score_result = WFAligner::builder(AlignmentScope::Score, MemoryModel::MemorySingletrack)
        .affine(6, 4, 2)
        .build();
    assert!(matches!(
        score_result,
        Err(WfaError::IncompatibleMemoryModel { .. })
    ));

    for result in [
        WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemorySingletrack)
            .indel()
            .build(),
        WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemorySingletrack)
            .edit()
            .build(),
        WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemorySingletrack)
            .linear(6, 2)
            .build(),
        WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemorySingletrack)
            .affine(6, 4, 2)
            .with_heuristics(Heuristics::banded_static(-10, 10))
            .build(),
        WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemorySingletrack)
            .affine(6, 4, 2)
            .with_heuristics(Heuristics::banded_adaptive(1, -10, 10))
            .build(),
    ] {
        assert!(matches!(
            result,
            Err(WfaError::IncompatibleMemoryModel { .. })
        ));
    }
}

#[test]
fn test_singletrack_set_heuristics_rejects_banded() {
    let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemorySingletrack)
        .affine(6, 4, 2)
        .build()
        .unwrap();

    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        aligner.set_heuristics(Heuristics::banded_static(-10, 10));
    }));

    let payload = result.expect_err("expected MemorySingletrack banded heuristic rejection");
    assert_eq!(
        panic_message(payload.as_ref()),
        Some("MemorySingletrack is incompatible with this aligner: singletrack does not support banded heuristics")
    );
}

#[test]
fn test_singletrack_rejects_lambda_inputs_before_entering_c() {
    let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemorySingletrack)
        .affine(6, 4, 2)
        .build()
        .unwrap();

    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        aligner.align_end_to_end_lambda(4, 4, |_, _| true);
    }));

    let payload = result.expect_err("expected MemorySingletrack lambda rejection");
    assert_eq!(
        panic_message(payload.as_ref()),
        Some("Lambda/custom sequence inputs are not supported with MemorySingletrack")
    );
}

#[test]
fn test_set_heuristics_replaces_configuration() {
    let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .affine(6, 4, 2)
        .build()
        .unwrap();
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
        .build()
        .unwrap();
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
        .build()
        .unwrap();
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
        .build()
        .unwrap();

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
        .build()
        .unwrap();

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
        .build()
        .unwrap();

    let result = aligner.align_end_to_end(PATTERN, TEXT);
    assert_eq!(result.status, AlignmentStatus::StatusMaxStepsReached);
    aligner.get_alignment();
}

#[test]
#[should_panic(expected = "Cannot get alignment when AlignmentScope is Score")]
fn test_get_alignment_rejects_score_scope() {
    let aligner = WFAligner::builder(AlignmentScope::Score, MemoryModel::MemoryHigh)
        .edit()
        .build()
        .unwrap();

    aligner.get_alignment();
}

#[test]
#[should_panic(expected = "Cannot get alignment span when AlignmentScope is Score")]
fn test_get_alignment_span_rejects_score_scope() {
    let aligner = WFAligner::builder(AlignmentScope::Score, MemoryModel::MemoryHigh)
        .edit()
        .build()
        .unwrap();

    aligner.get_alignment_span();
}

#[test]
fn ultralow_memory_default_heuristic_can_be_unattainable() {
    let read = b"GCTGCTACTGGGGTGTCCCCTCTCAAAGGACAAACCCAGGATCTACAGATGTGTGTGCTAAGCCATGTATGCACATGCACGTGTGTGTGTATATATTTAACCTATCTGTATATATGTATTATGTAAACATGAGTTCCTGCTGGCATATCTGACTATAACTGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACATGAGTTCCCTGCTGGCATATCTGACTATAACTGACCACCTCACAGTCCATTCTGATCTCTATATATGTATTATGTAAACATGAGTTCCTACTGGCATATCTGACTATAACTGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATTATGTAAACATGAGTTCCCTGCTGGCATATCTGATTATAACTGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATTATGTAAACATGAGTTCCTACTGGCATATCTGACTATAACCGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACACGAGTTCCTACTGGCATATCTGACTATAACTGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACACGAGTTCCTGCTGGCATATCTGACTATAACTGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATAATATATATTATATATGGACCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACATGAGTTCCTGCTGGCATATCTGACTATAACTGACCACCTCAGGGTCTATTCTGATCTGTATATATGTATAATATATATTATATATGGACCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACATGAGTTCCTGCTGGCATATCTGATTATAACCGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATTATGTAAACATGAGTTCCTACTGGCATATCTGACTATAACCGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACATGAGTTCCTACTGGCATATCTGACTATAACTGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACATGAGTTCCTACTGGCATATCTGACTATAACTGACCACCTCAGGATCCATTCTGATCTGTATATATGTATAATATATATTATATATGGACCTCAGGGTCCATTCTGATCTGTATATATGTATAATATATATTATATATGGACCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACATGAGTTCCTGGCTGGCATATCTGATTATAACCGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACACGAGTTCCTACTGGCATATCTGACTATAACTGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACACGAGTTCCTGCTGGCATATCTGATTATAACCGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATAATATATATTATATATGGACCTCAGGGTCCCCGCTGGCTTTTCCATGACTTCCTTATCCAGCTGTGAGAACCCTGACTCTTACTACCCATACTGTATTGACTTATTT";
    let allele = b"GCTGCTACTGGGGTGTCCCCTCTCAAAGGACAAACCCAGGATCTACAGATGTGTGTGCTAAGCCATGTATGCACACGCACGTGTGTGTGTATATATTTAACCTATCTGTATATATGTATTATGTAAACATGAGTTCCTGCTGGCATATCTGACTATAACTGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACACGACTTCCTACTGGCATATCTGACTGTAACCGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACATGATTTCCTACTGGCATATCTGACTATAACTGACCACCTCAGGGTTCATTCCGATCTGTATATAAGTATCATGTAAACACGAGTTCCTGCTGGCATATCTGACTGTAACCGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACACGAGTTCCTGCTGGCATATCTGACTATAACCGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACATGAGTTCCTACTGGCATATCTGACTATAACTGACCACCTCAGGGTCCATTCTGATCTGTATGTATGTATCATGTAAACACGAGTTCCTACTGGCATATCTGACTATAACTGACCACCTCAGGGTCCATTCCGATCTGTATATAAGTATCATGTAAACACGAGTTCCTGCTGGCATATCTGACTGTAACCGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACACGAGTTCCTGCTGGCATATCTGACTATAACTGACCACCTCAGGGTCCATTCTGATCTGTATATATGTATAATATATATTATATATGGACCTCAGGGTCCATTCTGATCTGCATATATGTATAATATATATTATATATGGACCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACATGAGTTCCTGCTGGCATATCTGTCTATAACCGACCACCTTAGGGTCCATTCTGATCTGTATATATGTATAATATATATTATATATGGTCCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACATGAGTTCCTGCTGGCATATCTGTCTATAACCGACCACCTTAGGGTCCATTCTGATCTGTATATATGTATAATATATATTATATATGGACCTCAGGGTCCATTCTGATCTGCATATATGTATAATATATATTATATATGGTCCTCAGGGTCCATTCTGATCTGTATATATGTATCATGTAAACATGAGTTCCTGCTGGCATATCTGTCTATAACCGACCACCTTAGGGTCCATTCTGATCTGTATATATGTATAATATATATTATATATGGACCTCAGGGTCCCCGCTGGCTTTTCCATGACTTCCTTATCCAGCTGTGAGAACCCTGACTCTTACTACTGTATTGACTTATTTGTGAAACCT";

    let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryUltraLow)
        .with_heuristics(Heuristics::wfa2_default())
        .affine2p(8, 4, 2, 24, 1)
        .build()
        .unwrap();

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
        .build()
        .unwrap();
    assert_eq!(aligner_edit.get_penalties(), Penalties::Edit);

    let aligner_indel = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
        .indel()
        .build()
        .unwrap();
    assert_eq!(aligner_indel.get_penalties(), Penalties::Indel);

    let aligner_linear = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
        .linear(12, 24)
        .build()
        .unwrap();
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
        .build()
        .unwrap();
    assert_eq!(
        aligner_affine.get_penalties(),
        Penalties::Affine {
            match_: 0,
            mismatch: 12,
            gap_opening: 24,
            gap_extension: 2
        }
    );

    let aligner_affine2p = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
        .affine2p(12, 24, 2, 48, 1)
        .build()
        .unwrap();
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
            .build()
            .unwrap();
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
        .build()
        .unwrap();
    assert_eq!(aligner_edit.get_penalties(), Penalties::Edit);

    let aligner_affine = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
        .affine(12, 24, 2)
        .with_heuristics(Heuristics::wf_adaptive(100, 10, 50))
        .build()
        .unwrap();
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
            .build()
            .unwrap();
    assert_eq!(
        aligner_affine_heuristic.get_penalties(),
        Penalties::Affine {
            match_: 0,
            mismatch: 12,
            gap_opening: 24,
            gap_extension: 2
        }
    );

    let aligner_affine2p = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
        .affine2p(12, 24, 2, 48, 1)
        .build()
        .unwrap();
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
        .build()
        .unwrap();
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
fn test_get_and_decode_packed_cigar() {
    let pattern = b"TCTTTACTCTT";
    let text = b"TCTTTACTCTT";
    let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .affine(4, 6, 2)
        .build()
        .unwrap();

    let result = aligner.align_end_to_end(pattern, text);
    assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);

    let sam_cigar_buffer = aligner.sam_packed_cigar(true);
    assert!(
        !sam_cigar_buffer.is_empty(),
        "SAM CIGAR buffer should not be empty"
    );

    let decoded_cigar = WFAligner::decode_packed_cigar(&sam_cigar_buffer);

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
    let sam_cigar_buffer_m = aligner.sam_packed_cigar(false);
    // 'M' is op_code 0. (11 << 4) | 0 = 176
    let expected_raw_buffer_m = vec![176]; // 11M
    assert_eq!(
        sam_cigar_buffer_m, expected_raw_buffer_m,
        "Raw SAM CIGAR buffer mismatch (M)"
    );

    let decoded_cigar_m = WFAligner::decode_packed_cigar(&sam_cigar_buffer_m);
    let expected_decoded_cigar_m: Vec<CigarOp> = vec![(11, 'M')]; // 11 matches ('M')
    assert_eq!(
        decoded_cigar_m, expected_decoded_cigar_m,
        "Decoded SAM CIGAR mismatch (M)"
    );

    let pattern_diff = b"TCTTTACTCTT";
    let text_diff = b"TCTTTACTATT";
    let mut aligner_diff = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
        .affine(4, 6, 2)
        .build()
        .unwrap();
    let result_diff = aligner_diff.align_end_to_end(pattern_diff, text_diff);
    assert_eq!(result_diff.status, AlignmentStatus::StatusAlgCompleted);

    let sam_cigar_buffer_diff = aligner_diff.sam_packed_cigar(true);

    let expected_raw_diff = vec![135, 24, 39];
    assert_eq!(
        sam_cigar_buffer_diff, expected_raw_diff,
        "Raw SAM CIGAR buffer mismatch (diff)"
    );

    let decoded_cigar_diff = WFAligner::decode_packed_cigar(&sam_cigar_buffer_diff);
    let expected_decoded_diff: Vec<CigarOp> = vec![(8, '='), (1, 'X'), (2, '=')];
    assert_eq!(
        decoded_cigar_diff, expected_decoded_diff,
        "Decoded SAM CIGAR mismatch (diff)"
    );

    // Test with show_mismatches = false
    let sam_cigar_buffer_diff_m = aligner_diff.sam_packed_cigar(false);
    // Expected: 11M => (11<<4)|0 = 176
    let expected_raw_diff_m = vec![176];
    assert_eq!(
        sam_cigar_buffer_diff_m, expected_raw_diff_m,
        "Raw SAM CIGAR buffer mismatch (diff, M)"
    );

    let decoded_cigar_diff_m = WFAligner::decode_packed_cigar(&sam_cigar_buffer_diff_m);
    let expected_decoded_diff_m: Vec<CigarOp> = vec![(11, 'M')];
    assert_eq!(
        decoded_cigar_diff_m, expected_decoded_diff_m,
        "Decoded SAM CIGAR mismatch (diff, M)"
    );
}

#[test]
fn test_wfa_and_sam_cigars_have_explicit_indel_orientation() {
    let query = b"ACGTT";
    let reference = b"ACGT";
    let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .edit()
        .build()
        .unwrap();

    let result = aligner.align_end_to_end(query, reference);
    assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
    assert_eq!(aligner.wfa_cigar_bytes(), b"MMMMD");
    assert_eq!(aligner.sam_cigar_bytes(), b"MMMMI");
    assert_eq!(aligner.wfa_packed_cigar(false), vec![64, 18]);
    assert_eq!(aligner.sam_packed_cigar(false), vec![64, 17]);
    assert_eq!(aligner.wfa_cigar(false), vec![(4, 'M'), (1, 'D')]);
    assert_eq!(aligner.sam_cigar(false), vec![(4, 'M'), (1, 'I')]);

    let result = aligner.align_end_to_end(reference, query);
    assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
    assert_eq!(aligner.wfa_cigar_bytes(), b"MMMMI");
    assert_eq!(aligner.wfa_packed_cigar(false), vec![64, 17]);
}

#[test]
fn test_get_heuristics_round_trips_combined_categories() {
    let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .affine(12, 24, 2)
        .build()
        .unwrap();
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
            .build()
            .unwrap();
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
        .build()
        .unwrap();

    aligner.set_heuristics(Heuristics::none());
    assert!(aligner.get_heuristics().is_none());
    assert_eq!(aligner.get_heuristics(), Heuristics::none());
}

#[test]
fn test_drop_heuristics_reject_edit_and_indel() {
    let edit_result = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .edit()
        .with_heuristics(Heuristics::xdrop(1, 10))
        .build();
    assert!(matches!(
        edit_result,
        Err(WfaError::IncompatibleHeuristics { .. })
    ));

    let indel_result = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .indel()
        .with_heuristics(Heuristics::zdrop(1, 10))
        .build();
    assert!(matches!(
        indel_result,
        Err(WfaError::IncompatibleHeuristics { .. })
    ));

    let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .edit()
        .build()
        .unwrap();
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
        .build()
        .unwrap();

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
        .build()
        .unwrap();

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
fn test_cigar_view_active_operations_reject_invalid_offsets() {
    let operations = b"MMIDX"
        .iter()
        .map(|&operation| operation as c_char)
        .collect::<Vec<_>>();

    let valid_cigar = CigarView::new(0, 1, 4, 0, 0, &operations);
    assert_eq!(valid_cigar.active_operation_bytes(), b"MID");

    let invalid_ranges = [(-1, 3), (0, 6), (4, 3)];

    for (begin_offset, end_offset) in invalid_ranges {
        let cigar = CigarView::new(0, begin_offset, end_offset, 0, 0, &operations);

        assert!(cigar.active_operation_bytes().is_empty());
    }
}

#[test]
fn test_get_alignment_global() {
    let pattern = b"AGCTAGTGTCAATGGCTACTTTTCAGGTCCT";
    let text = b"AACTAAGTGTCGGTGGCTACTATATATCAGGTCCT";
    let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .affine(1, 5, 1)
        .build()
        .unwrap();
    let result = aligner.align_end_to_end(pattern, text);
    assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);

    let alignment = aligner.get_alignment();
    assert_eq!(aligner.score(), -18);

    let expected_ops = vec![
        Match, Subst, Match, Match, Match, Ins, Match, Match, Match, Match, Match, Subst, Subst,
        Match, Match, Match, Match, Match, Match, Match, Match, Ins, Ins, Ins, Match, Subst, Match,
        Match, Match, Match, Match, Match, Match, Match, Match,
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
    let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryUltraLow)
        .affine(1, 5, 1)
        .build()
        .unwrap();
    let result = aligner.align_end_to_end(pattern, text);
    assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);

    let alignment = aligner.get_alignment();
    assert_eq!(aligner.score(), -18);
    assert_eq!(aligner.cigar_score(), -18);

    let expected_ops = vec![
        Match, Subst, Match, Match, Match, Ins, Match, Match, Match, Match, Match, Subst, Subst,
        Match, Match, Match, Match, Match, Match, Match, Match, Ins, Ins, Ins, Match, Subst, Match,
        Match, Match, Match, Match, Match, Match, Match, Match,
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

    let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryUltraLow)
        .affine(1, 5, 1)
        .build()
        .unwrap();
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
        .build()
        .unwrap();
    let result = aligner.align_ends_free(pattern, 0, 0, text, text.len() as i32, text.len() as i32);
    assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
    let alignment = aligner.get_alignment();
    assert_eq!(aligner.score(), 0);

    let expected_ops = vec![
        Ins, Ins, Ins, Ins, Ins, Ins, Ins, Ins, Ins, Ins, Match, Match, Match, Match, Match, Match,
        Match, Match, Match, Match, Match, Match, Match, Match, Match, Ins, Ins, Ins, Ins, Ins,
        Ins, Ins, Ins, Ins, Ins,
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

#[test]
#[should_panic(expected = "max_alignment_steps must be positive")]
fn test_resource_limits_rejects_nonpositive_max_alignment_steps() {
    ResourceLimits::new(0, 100, 100, 1, 1);
}

#[test]
#[should_panic(expected = "max_memory_resident must be less than or equal to max_memory_abort")]
fn test_resource_limits_rejects_resident_above_abort() {
    ResourceLimits::new(1, 100, 50, 1, 1);
}

#[test]
#[should_panic(expected = "max_num_threads must be positive")]
fn test_resource_limits_rejects_nonpositive_threads() {
    ResourceLimits::new(1, 50, 50, 0, 1);
}

#[test]
#[should_panic(expected = "min_offsets_per_thread must be positive")]
fn test_resource_limits_rejects_nonpositive_min_offsets() {
    ResourceLimits::new(1, 50, 50, 1, 0);
}

#[test]
fn test_resource_limits_allows_equal_memory_thresholds() {
    // The resident <= abort invariant is inclusive at the boundary.
    let limits = ResourceLimits::new(1, 50, 50, 1, 1);
    assert_eq!(limits.max_memory_resident, 50);
    assert_eq!(limits.max_memory_abort, 50);
}

#[test]
#[should_panic(expected = "max_memory_resident must be less than or equal to max_memory_abort")]
fn test_builder_with_max_memory_validates() {
    let _ = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .with_max_memory(100, 50);
}

#[test]
#[should_panic(expected = "max_alignment_steps must be positive")]
fn test_runtime_set_max_alignment_steps_validates() {
    let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .edit()
        .build()
        .unwrap();
    aligner.set_max_alignment_steps(0);
}

#[test]
#[should_panic(expected = "resolution_points must be positive")]
fn test_plot_options_rejects_nonpositive_resolution() {
    PlotOptions::new(0, 0);
}

#[test]
#[should_panic(expected = "align_level must be greater than or equal to -1")]
fn test_plot_options_rejects_align_level_below_minus_one() {
    PlotOptions::new(2000, -2);
}

#[test]
fn test_plot_options_allows_final_alignment_sentinel() {
    // align_level == -1 is the valid "final/subsidiary alignment" sentinel.
    let options = PlotOptions::new(2000, -1);
    assert_eq!(options.align_level, -1);
    assert_eq!(PlotOptions::final_alignment().align_level, -1);
    assert_eq!(PlotOptions::at_recursion_level(0).align_level, 0);
}

#[test]
#[should_panic(expected = "steps_between_cutoffs must be positive")]
fn test_heuristics_rejects_nonpositive_steps() {
    Heuristics::new(0);
}

#[test]
#[should_panic(expected = "min_wavefront_length must be positive")]
fn test_heuristics_rejects_nonpositive_min_wavefront_length() {
    Heuristics::wf_adaptive(1, 0, 50);
}

#[test]
#[should_panic(expected = "max_distance_threshold must be non-negative")]
fn test_heuristics_rejects_negative_max_distance_threshold() {
    Heuristics::wf_adaptive(1, 1, -1);
}

#[test]
fn test_heuristics_allows_zero_max_distance_threshold() {
    let heuristics = Heuristics::wf_adaptive(1, 1, 0);
    assert_eq!(
        heuristics.adaptive(),
        Some(AdaptiveHeuristic::WfAdaptive {
            min_wavefront_length: 1,
            max_distance_threshold: 0,
        })
    );
}

#[test]
#[should_panic(expected = "xdrop must be non-negative")]
fn test_heuristics_rejects_negative_xdrop() {
    Heuristics::xdrop(1, -1);
}

#[test]
#[should_panic(expected = "zdrop must be non-negative")]
fn test_heuristics_rejects_negative_zdrop() {
    Heuristics::zdrop(1, -1);
}

#[test]
fn test_heuristics_allows_zero_drop() {
    assert_eq!(
        Heuristics::xdrop(1, 0).drop_heuristic(),
        Some(DropHeuristic::XDrop { xdrop: 0 })
    );
}

#[test]
#[should_panic(expected = "min_k must be less than or equal to max_k")]
fn test_heuristics_rejects_inverted_band() {
    Heuristics::banded_static(5, 4);
}

#[test]
fn test_heuristics_allows_equal_band_bounds() {
    assert_eq!(
        Heuristics::banded_static(5, 5).band(),
        Some(BandHeuristic::Static { min_k: 5, max_k: 5 })
    );
}

#[test]
fn test_builder_rejects_missing_penalty_model() {
    match WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh).build() {
        Err(err) => assert_eq!(err, WfaError::MissingPenaltyModel),
        Ok(_) => panic!("expected missing penalty model error"),
    }
}

#[test]
#[allow(clippy::identity_op)] // `| 0` kept for column alignment with the other op codes
fn test_decode_packed_cigar_covers_all_op_codes_and_unknown_fallback() {
    // Encoding is (length << 4) | op_code. Cover every documented op code plus an
    // out-of-range code (>= 9) that must decode to the '?' fallback.
    let buffer = vec![
        (1u32 << 4) | 0,  // 1M
        (2u32 << 4) | 1,  // 2I
        (3u32 << 4) | 2,  // 3D
        (4u32 << 4) | 3,  // 4N
        (5u32 << 4) | 4,  // 5S
        (6u32 << 4) | 5,  // 6H
        (7u32 << 4) | 6,  // 7P
        (8u32 << 4) | 7,  // 8=
        (9u32 << 4) | 8,  // 9X
        (10u32 << 4) | 9, // 10? (unknown op code)
    ];
    let decoded = WFAligner::decode_packed_cigar(&buffer);
    assert_eq!(
        decoded,
        vec![
            (1, 'M'),
            (2, 'I'),
            (3, 'D'),
            (4, 'N'),
            (5, 'S'),
            (6, 'H'),
            (7, 'P'),
            (8, '='),
            (9, 'X'),
            (10, '?'),
        ]
    );
}

#[test]
fn test_decode_packed_cigar_empty_buffer() {
    assert!(WFAligner::decode_packed_cigar(&[]).is_empty());
}

#[test]
#[allow(clippy::identity_op)] // `| 0` documents the 'M' op code (0) explicitly
fn test_decode_packed_cigar_decodes_large_lengths() {
    // 28-bit maximum length must survive the >> 4 shift without truncation.
    let max_len = (1u32 << 28) - 1;
    let buffer = vec![(max_len << 4) | 0];
    assert_eq!(
        WFAligner::decode_packed_cigar(&buffer),
        vec![(max_len as usize, 'M')]
    );
}

#[test]
#[should_panic(expected = "Cannot clip when AlignmentScope is Score")]
fn test_cigar_score_clipped_rejects_score_scope() {
    let aligner = WFAligner::builder(AlignmentScope::Score, MemoryModel::MemoryHigh)
        .edit()
        .build()
        .unwrap();
    aligner.cigar_score_clipped(0);
}

#[test]
#[should_panic(expected = "Cannot count matches when AlignmentScope is Score")]
fn test_count_matches_rejects_score_scope() {
    let aligner = WFAligner::builder(AlignmentScope::Score, MemoryModel::MemoryHigh)
        .edit()
        .build()
        .unwrap();
    aligner.count_matches();
}

#[test]
#[should_panic(expected = "Cannot get SAM CIGAR when AlignmentScope is Score")]
fn test_sam_packed_cigar_rejects_score_scope() {
    let aligner = WFAligner::builder(AlignmentScope::Score, MemoryModel::MemoryHigh)
        .edit()
        .build()
        .unwrap();
    aligner.sam_packed_cigar(true);
}

#[test]
#[should_panic(expected = "Cannot get WFA packed CIGAR when AlignmentScope is Score")]
fn test_wfa_packed_cigar_rejects_score_scope() {
    let aligner = WFAligner::builder(AlignmentScope::Score, MemoryModel::MemoryHigh)
        .edit()
        .build()
        .unwrap();
    aligner.wfa_packed_cigar(true);
}

#[test]
#[should_panic(expected = "Cannot calculate CIGAR score when AlignmentScope is Score")]
fn test_cigar_score_rejects_score_scope() {
    let mut aligner = WFAligner::builder(AlignmentScope::Score, MemoryModel::MemoryHigh)
        .edit()
        .build()
        .unwrap();
    aligner.cigar_score();
}

#[test]
#[should_panic(expected = "Cannot get WFA CIGAR bytes when AlignmentScope is Score")]
fn test_wfa_cigar_bytes_rejects_score_scope() {
    let mut aligner = WFAligner::builder(AlignmentScope::Score, MemoryModel::MemoryHigh)
        .edit()
        .build()
        .unwrap();
    let result = aligner.align_end_to_end(PATTERN, TEXT);
    assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
    aligner.wfa_cigar_bytes();
}

#[test]
#[should_panic(expected = "Cannot get WFA CIGAR bytes when AlignmentScope is Score")]
fn test_sam_cigar_bytes_rejects_score_scope() {
    let mut aligner = WFAligner::builder(AlignmentScope::Score, MemoryModel::MemoryHigh)
        .edit()
        .build()
        .unwrap();
    let result = aligner.align_end_to_end(PATTERN, TEXT);
    assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
    aligner.sam_cigar_bytes();
}

#[test]
#[should_panic(expected = "Unknown alignment status")]
fn test_alignment_status_from_unknown_panics() {
    let _ = AlignmentStatus::from(123_456_i32);
}

#[test]
#[should_panic(expected = "Unknown distance metric")]
fn test_distance_metric_from_unknown_panics() {
    let _ = DistanceMetric::from(9999_u32);
}

#[test]
#[should_panic(expected = "Unknown alignment scope")]
fn test_alignment_scope_from_unknown_panics() {
    let _ = AlignmentScope::from(9999_u32);
}

#[test]
#[should_panic(expected = "Invalid alignment operation character")]
fn test_wfa_op_from_invalid_byte_panics() {
    let _ = WfaOp::from_u8(b'Z');
}

#[test]
fn test_alignment_status_i32_round_trip() {
    // Guards against drift between the enum discriminants and the From<i32> mapping.
    for status in [
        AlignmentStatus::StatusAlgCompleted,
        AlignmentStatus::StatusAlgPartial,
        AlignmentStatus::StatusMaxStepsReached,
        AlignmentStatus::StatusOOM,
        AlignmentStatus::StatusUnattainable,
    ] {
        assert_eq!(AlignmentStatus::from(status as i32), status);
    }
}

#[test]
fn test_alignment_status_display_strings() {
    assert_eq!(
        format!("{}", AlignmentStatus::StatusAlgCompleted),
        "StatusAlgCompleted"
    );
    assert_eq!(
        format!("{}", AlignmentStatus::StatusAlgPartial),
        "StatusAlgPartial"
    );
    assert_eq!(
        format!("{}", AlignmentStatus::StatusMaxStepsReached),
        "StatusMaxStepsReached"
    );
    assert_eq!(format!("{}", AlignmentStatus::StatusOOM), "StatusOOM");
    assert_eq!(
        format!("{}", AlignmentStatus::StatusUnattainable),
        "StatusUnattainable"
    );
}

fn divergent_sequences(len: usize) -> (Vec<u8>, Vec<u8>) {
    // Two independent pseudo-random sequences. Being unrelated, they force a high score and
    // therefore a large MemoryHigh wavefront footprint.
    let bases = [b'A', b'C', b'G', b'T'];
    let mut pattern = Vec::with_capacity(len);
    let mut text = Vec::with_capacity(len);
    let mut state: u64 = 0x1234_5678_9ABC_DEF0;
    let mut next_base = || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        bases[((state >> 33) % 4) as usize]
    };
    for _ in 0..len {
        pattern.push(next_base());
        text.push(next_base());
    }
    (pattern, text)
}

#[test]
fn test_alignment_aborts_with_oom_under_tiny_memory_budget() {
    let (pattern, text) = divergent_sequences(2000);
    let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .with_max_memory(1024, 1024)
        .affine(6, 4, 2)
        .build()
        .unwrap();

    let result = aligner.align_end_to_end(&pattern, &text);
    assert_eq!(result.status, AlignmentStatus::StatusOOM);
}

#[test]
#[should_panic]
fn test_get_alignment_span_rejects_oom() {
    let (pattern, text) = divergent_sequences(2000);
    let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .with_max_memory(1024, 1024)
        .affine(6, 4, 2)
        .build()
        .unwrap();

    let result = aligner.align_end_to_end(&pattern, &text);
    assert_eq!(result.status, AlignmentStatus::StatusOOM);
    aligner.get_alignment_span();
}

#[test]
fn test_write_plot_rejects_interior_nul_in_path() {
    let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .with_plotting(PlotOptions::default())
        .edit()
        .build()
        .unwrap();
    let result = aligner.align_end_to_end(PATTERN, TEXT);
    assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);

    let err = aligner.write_plot(Path::new("bad\0name.plot")).unwrap_err();
    assert_eq!(err.kind(), io::ErrorKind::InvalidInput);
}

#[test]
fn test_align_empty_pattern_against_text() {
    let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .affine(6, 4, 2)
        .build()
        .unwrap();

    let text = b"ACGT";
    let result = aligner.align_end_to_end(b"", text);
    assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
    assert_eq!(raw_cigar_string(&aligner), "IIII");

    let alignment = aligner.get_alignment();
    assert_eq!(alignment.xlen, 0);
    assert_eq!(alignment.ylen, text.len());
    assert_eq!(alignment.operations, vec![Ins, Ins, Ins, Ins]);
    assert_eq!(aligner.get_alignment_span(), ((0, 0), (0, 4)));
}

#[test]
fn test_align_both_empty() {
    let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .affine(6, 4, 2)
        .build()
        .unwrap();

    let result = aligner.align_end_to_end(b"", b"");
    assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
    assert_eq!(aligner.score(), 0);

    let alignment = aligner.get_alignment();
    assert!(alignment.operations.is_empty());
    assert_eq!(alignment.xlen, 0);
    assert_eq!(alignment.ylen, 0);
    assert_eq!(aligner.get_alignment_span(), ((0, 0), (0, 0)));
    assert!(aligner.wfa_cigar_bytes().is_empty());
}

#[test]
fn test_ultralow_ends_free_all_zero_matches_global() {
    let pattern = b"AGCTAGTGTCAATGGCTACTTTTCAGGTCCT";
    let text = b"AACTAAGTGTCGGTGGCTACTATATATCAGGTCCT";

    let mut global = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryUltraLow)
        .affine(1, 5, 1)
        .build()
        .unwrap();
    let global_result = global.align_end_to_end(pattern, text);
    assert_eq!(global_result.status, AlignmentStatus::StatusAlgCompleted);
    let global_score = global.score();

    // All-zero free ends are permitted with MemoryUltraLow (the guard only fires on a
    // nonzero free end) and must degenerate to the global alignment.
    let mut ends_free = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryUltraLow)
        .affine(1, 5, 1)
        .build()
        .unwrap();
    let ends_free_result = ends_free.align_ends_free(pattern, 0, 0, text, 0, 0);
    assert_eq!(ends_free_result.status, AlignmentStatus::StatusAlgCompleted);
    assert_eq!(ends_free.score(), global_score);
}

#[test]
fn test_cigar_score_clipped_flank_exceeding_alignment_is_zero() {
    let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
        .affine2p_with_match(-1, 3, 3, 3, 10, 0)
        .build()
        .unwrap();

    let pattern = b"TCTATAATAGT";
    let text = b"TCTATACTGCGCGTTTGGAGAAATAAAATAGT";
    let result = aligner.align_end_to_end(pattern, text);
    assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
    assert_eq!(
        raw_cigar_string(&aligner),
        "MMMMMMIIIIIIIIIIIIIIIIIIIIIMMMMM"
    );

    // A flank that meets or exceeds half the CIGAR collapses the clipped window to empty.
    assert_eq!(aligner.cigar_score_clipped(1000), 0);
    assert_eq!(aligner.cigar_score_clipped(16), 0);
}

#[test]
fn test_cigar_score_clipped_affine2p_selects_cheaper_gap_piece() {
    let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
        .affine2p_with_match(-1, 3, 3, 3, 10, 0)
        .build()
        .unwrap();

    let pattern = b"TCTATAATAGT";
    let text = b"TCTATACTGCGCGTTTGGAGAAATAAAATAGT";
    let result = aligner.align_end_to_end(pattern, text);
    assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);

    // CIGAR is 6M21I5M. The 21I gap must be scored by the cheaper second piece
    // (10 + 0*21 = 10), not the first (3 + 3*21 = 66): 6*(-1) negated = 6 for the matches,
    // minus 10 for the gap -> 1. Selecting piece 1 would give 6 - 66 = -55.
    assert_eq!(aligner.cigar_score_clipped(0), 1);

    // Clipping into the middle leaves a pure 16I window, also scored by piece 2
    // (10 + 0*16 = 10).
    assert_eq!(aligner.cigar_score_clipped(8), -10);
}

#[test]
fn test_count_matches_direct() {
    let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .affine(4, 6, 2)
        .build()
        .unwrap();

    let identical = b"TCTTTACTCTT";
    let result = aligner.align_end_to_end(identical, identical);
    assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
    assert_eq!(aligner.count_matches(), 11);

    let with_mismatch = b"TCTTTACTATT";
    let result = aligner.align_end_to_end(identical, with_mismatch);
    assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
    assert_eq!(aligner.count_matches(), 10);
}

#[test]
fn test_heuristic_constructors_drive_completed_alignments() {
    for heuristics in [
        Heuristics::wfa2_default(),
        Heuristics::wf_mash(1, 10, 50),
        Heuristics::banded_adaptive(1, -50, 50),
    ] {
        let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
            .affine(6, 4, 2)
            .with_heuristics(heuristics)
            .build()
            .unwrap();
        let result = aligner.align_end_to_end(PATTERN, TEXT);
        assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
    }
}

#[cfg(feature = "serde")]
#[test]
fn test_serde_round_trips_config_types() {
    let penalties = Penalties::Affine2p {
        match_: -1,
        mismatch: 8,
        gap_opening1: 4,
        gap_extension1: 2,
        gap_opening2: 24,
        gap_extension2: 1,
    };
    let json = serde_json::to_string(&penalties).unwrap();
    assert_eq!(serde_json::from_str::<Penalties>(&json).unwrap(), penalties);

    let heuristics = Heuristics::new(3)
        .with_adaptive(AdaptiveHeuristic::WfMash {
            min_wavefront_length: 5,
            max_distance_threshold: 25,
        })
        .with_drop(DropHeuristic::ZDrop { zdrop: 15 })
        .with_band(BandHeuristic::Static {
            min_k: -10,
            max_k: 10,
        });
    let json = serde_json::to_string(&heuristics).unwrap();
    assert_eq!(
        serde_json::from_str::<Heuristics>(&json).unwrap(),
        heuristics
    );

    let limits = ResourceLimits::new(64, 1_048_576, 2_097_152, 1, 64);
    let json = serde_json::to_string(&limits).unwrap();
    assert_eq!(
        serde_json::from_str::<ResourceLimits>(&json).unwrap(),
        limits
    );

    let plot = PlotOptions::new(1500, -1);
    let json = serde_json::to_string(&plot).unwrap();
    assert_eq!(serde_json::from_str::<PlotOptions>(&json).unwrap(), plot);
}

#[cfg(feature = "serde")]
#[test]
fn test_serde_round_trips_result_types() {
    let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .affine(6, 4, 2)
        .build()
        .unwrap();
    let result = aligner.align_end_to_end(PATTERN, TEXT);
    assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);

    let json = serde_json::to_string(&result).unwrap();
    assert_eq!(
        serde_json::from_str::<AlignmentResult>(&json).unwrap(),
        result
    );

    let alignment = aligner.get_alignment();
    let json = serde_json::to_string(&alignment).unwrap();
    assert_eq!(serde_json::from_str::<WfaAlign>(&json).unwrap(), alignment);
}

#[cfg(feature = "serde")]
#[test]
#[should_panic(expected = "min_k must be less than or equal to max_k")]
fn test_serde_deserialized_invalid_config_is_revalidated_on_use() {
    // Deserialization itself does not run the constructor validators, but every config
    // type is re-validated at the FFI boundary, so an invalid deserialized value still
    // cannot reach WFA2 silently.
    let json = r#"{"steps_between_cutoffs":1,"adaptive":null,"drop_heuristic":null,"band":{"Static":{"min_k":10,"max_k":-10}}}"#;
    let heuristics: Heuristics = serde_json::from_str(json).unwrap();
    let _ = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
        .affine(6, 4, 2)
        .with_heuristics(heuristics);
}
