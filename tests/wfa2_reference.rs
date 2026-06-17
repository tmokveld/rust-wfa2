use rust_wfa2::aligner::{AlignmentScope, AlignmentStatus, Heuristics, MemoryModel, WFAligner};

const SEQUENCES: &str = include_str!("../wfa2-sys/WFA2-lib/tests/wfa.utest.seq");

#[derive(Clone, Copy)]
struct Case {
    name: &'static str,
    scope: AlignmentScope,
    memory: MemoryModel,
    penalties: Penalties,
    heuristic: Option<ReferenceHeuristic>,
    expected: &'static str,
}

#[derive(Clone, Copy)]
enum Penalties {
    Indel,
    Edit,
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

#[derive(Clone, Copy)]
enum ReferenceHeuristic {
    WfAdaptive {
        min_wavefront_length: i32,
        max_distance_threshold: i32,
        steps_between_cutoffs: i32,
    },
}

#[derive(Debug)]
struct ExpectedRow<'a> {
    score: i32,
    cigar: &'a str,
}

macro_rules! reference_case {
    ($name:literal, $scope:expr, $memory:expr, $penalties:expr, $heuristic:expr) => {
        Case {
            name: $name,
            scope: $scope,
            memory: $memory,
            penalties: $penalties,
            heuristic: $heuristic,
            expected: include_str!(concat!(
                "../wfa2-sys/WFA2-lib/tests/wfa.utest.check/",
                $name,
                ".alg"
            )),
        }
    };
}

const AFFINE_DEFAULT: Penalties = Penalties::Affine {
    match_: 0,
    mismatch: 4,
    gap_opening: 6,
    gap_extension: 2,
};
const AFFINE2P_DEFAULT: Penalties = Penalties::Affine2p {
    match_: 0,
    mismatch: 4,
    gap_opening1: 6,
    gap_extension1: 2,
    gap_opening2: 24,
    gap_extension2: 1,
};
const AFFINE_P0: Penalties = Penalties::Affine {
    match_: 0,
    mismatch: 1,
    gap_opening: 2,
    gap_extension: 1,
};
const AFFINE_P1: Penalties = Penalties::Affine {
    match_: 0,
    mismatch: 3,
    gap_opening: 1,
    gap_extension: 4,
};
const AFFINE_P2: Penalties = Penalties::Affine {
    match_: 0,
    mismatch: 5,
    gap_opening: 3,
    gap_extension: 2,
};
const AFFINE_P3: Penalties = Penalties::Affine {
    match_: -5,
    mismatch: 1,
    gap_opening: 2,
    gap_extension: 1,
};
const AFFINE_P4: Penalties = Penalties::Affine {
    match_: -2,
    mismatch: 3,
    gap_opening: 1,
    gap_extension: 4,
};
const AFFINE_P5: Penalties = Penalties::Affine {
    match_: -3,
    mismatch: 5,
    gap_opening: 3,
    gap_extension: 2,
};
const WFAPT0: Option<ReferenceHeuristic> = Some(ReferenceHeuristic::WfAdaptive {
    min_wavefront_length: 10,
    max_distance_threshold: 50,
    steps_between_cutoffs: 1,
});
const WFAPT1: Option<ReferenceHeuristic> = Some(ReferenceHeuristic::WfAdaptive {
    min_wavefront_length: 10,
    max_distance_threshold: 50,
    steps_between_cutoffs: 10,
});

const CASES: &[Case] = &[
    reference_case!(
        "test.indel",
        AlignmentScope::Alignment,
        MemoryModel::MemoryHigh,
        Penalties::Indel,
        None
    ),
    reference_case!(
        "test.edit",
        AlignmentScope::Alignment,
        MemoryModel::MemoryHigh,
        Penalties::Edit,
        None
    ),
    reference_case!(
        "test.affine",
        AlignmentScope::Alignment,
        MemoryModel::MemoryHigh,
        AFFINE_DEFAULT,
        None
    ),
    reference_case!(
        "test.affine2p",
        AlignmentScope::Alignment,
        MemoryModel::MemoryHigh,
        AFFINE2P_DEFAULT,
        None
    ),
    reference_case!(
        "test.affine.p0",
        AlignmentScope::Alignment,
        MemoryModel::MemoryHigh,
        AFFINE_P0,
        None
    ),
    reference_case!(
        "test.affine.p1",
        AlignmentScope::Alignment,
        MemoryModel::MemoryHigh,
        AFFINE_P1,
        None
    ),
    reference_case!(
        "test.affine.p2",
        AlignmentScope::Alignment,
        MemoryModel::MemoryHigh,
        AFFINE_P2,
        None
    ),
    reference_case!(
        "test.affine.p3",
        AlignmentScope::Alignment,
        MemoryModel::MemoryHigh,
        AFFINE_P3,
        None
    ),
    reference_case!(
        "test.affine.p4",
        AlignmentScope::Alignment,
        MemoryModel::MemoryHigh,
        AFFINE_P4,
        None
    ),
    reference_case!(
        "test.affine.p5",
        AlignmentScope::Alignment,
        MemoryModel::MemoryHigh,
        AFFINE_P5,
        None
    ),
    reference_case!(
        "test.affine.wfapt0",
        AlignmentScope::Alignment,
        MemoryModel::MemoryHigh,
        AFFINE_DEFAULT,
        WFAPT0
    ),
    reference_case!(
        "test.affine.wfapt1",
        AlignmentScope::Alignment,
        MemoryModel::MemoryHigh,
        AFFINE_DEFAULT,
        WFAPT1
    ),
    reference_case!(
        "test.score.indel",
        AlignmentScope::Score,
        MemoryModel::MemoryHigh,
        Penalties::Indel,
        None
    ),
    reference_case!(
        "test.score.edit",
        AlignmentScope::Score,
        MemoryModel::MemoryHigh,
        Penalties::Edit,
        None
    ),
    reference_case!(
        "test.score.affine",
        AlignmentScope::Score,
        MemoryModel::MemoryHigh,
        AFFINE_DEFAULT,
        None
    ),
    reference_case!(
        "test.score.affine2p",
        AlignmentScope::Score,
        MemoryModel::MemoryHigh,
        AFFINE2P_DEFAULT,
        None
    ),
    reference_case!(
        "test.score.affine.p0",
        AlignmentScope::Score,
        MemoryModel::MemoryHigh,
        AFFINE_P0,
        None
    ),
    reference_case!(
        "test.score.affine.p1",
        AlignmentScope::Score,
        MemoryModel::MemoryHigh,
        AFFINE_P1,
        None
    ),
    reference_case!(
        "test.score.affine.p2",
        AlignmentScope::Score,
        MemoryModel::MemoryHigh,
        AFFINE_P2,
        None
    ),
    reference_case!(
        "test.score.affine.p3",
        AlignmentScope::Score,
        MemoryModel::MemoryHigh,
        AFFINE_P3,
        None
    ),
    reference_case!(
        "test.score.affine.p4",
        AlignmentScope::Score,
        MemoryModel::MemoryHigh,
        AFFINE_P4,
        None
    ),
    reference_case!(
        "test.score.affine.p5",
        AlignmentScope::Score,
        MemoryModel::MemoryHigh,
        AFFINE_P5,
        None
    ),
    reference_case!(
        "test.score.affine.wfapt0",
        AlignmentScope::Score,
        MemoryModel::MemoryHigh,
        AFFINE_DEFAULT,
        WFAPT0
    ),
    reference_case!(
        "test.score.affine.wfapt1",
        AlignmentScope::Score,
        MemoryModel::MemoryHigh,
        AFFINE_DEFAULT,
        WFAPT1
    ),
    reference_case!(
        "test.pb.indel",
        AlignmentScope::Alignment,
        MemoryModel::MemoryMed,
        Penalties::Indel,
        None
    ),
    reference_case!(
        "test.pb.edit",
        AlignmentScope::Alignment,
        MemoryModel::MemoryMed,
        Penalties::Edit,
        None
    ),
    reference_case!(
        "test.pb.affine",
        AlignmentScope::Alignment,
        MemoryModel::MemoryMed,
        AFFINE_DEFAULT,
        None
    ),
    reference_case!(
        "test.pb.affine2p",
        AlignmentScope::Alignment,
        MemoryModel::MemoryMed,
        AFFINE2P_DEFAULT,
        None
    ),
    reference_case!(
        "test.pb.affine.p0",
        AlignmentScope::Alignment,
        MemoryModel::MemoryMed,
        AFFINE_P0,
        None
    ),
    reference_case!(
        "test.pb.affine.p1",
        AlignmentScope::Alignment,
        MemoryModel::MemoryMed,
        AFFINE_P1,
        None
    ),
    reference_case!(
        "test.pb.affine.p2",
        AlignmentScope::Alignment,
        MemoryModel::MemoryMed,
        AFFINE_P2,
        None
    ),
    reference_case!(
        "test.pb.affine.p3",
        AlignmentScope::Alignment,
        MemoryModel::MemoryMed,
        AFFINE_P3,
        None
    ),
    reference_case!(
        "test.pb.affine.p4",
        AlignmentScope::Alignment,
        MemoryModel::MemoryMed,
        AFFINE_P4,
        None
    ),
    reference_case!(
        "test.pb.affine.p5",
        AlignmentScope::Alignment,
        MemoryModel::MemoryMed,
        AFFINE_P5,
        None
    ),
    reference_case!(
        "test.pb.affine.wfapt0",
        AlignmentScope::Alignment,
        MemoryModel::MemoryMed,
        AFFINE_DEFAULT,
        WFAPT0
    ),
    reference_case!(
        "test.pb.affine.wfapt1",
        AlignmentScope::Alignment,
        MemoryModel::MemoryMed,
        AFFINE_DEFAULT,
        WFAPT1
    ),
    reference_case!(
        "test.biwfa.indel",
        AlignmentScope::Alignment,
        MemoryModel::MemoryUltraLow,
        Penalties::Indel,
        None
    ),
    reference_case!(
        "test.biwfa.edit",
        AlignmentScope::Alignment,
        MemoryModel::MemoryUltraLow,
        Penalties::Edit,
        None
    ),
    reference_case!(
        "test.biwfa.affine",
        AlignmentScope::Alignment,
        MemoryModel::MemoryUltraLow,
        AFFINE_DEFAULT,
        None
    ),
    reference_case!(
        "test.biwfa.affine2p",
        AlignmentScope::Alignment,
        MemoryModel::MemoryUltraLow,
        AFFINE2P_DEFAULT,
        None
    ),
    reference_case!(
        "test.biwfa.affine.p0",
        AlignmentScope::Alignment,
        MemoryModel::MemoryUltraLow,
        AFFINE_P0,
        None
    ),
    reference_case!(
        "test.biwfa.affine.p1",
        AlignmentScope::Alignment,
        MemoryModel::MemoryUltraLow,
        AFFINE_P1,
        None
    ),
    reference_case!(
        "test.biwfa.affine.p2",
        AlignmentScope::Alignment,
        MemoryModel::MemoryUltraLow,
        AFFINE_P2,
        None
    ),
    reference_case!(
        "test.biwfa.affine.p3",
        AlignmentScope::Alignment,
        MemoryModel::MemoryUltraLow,
        AFFINE_P3,
        None
    ),
    reference_case!(
        "test.biwfa.affine.p4",
        AlignmentScope::Alignment,
        MemoryModel::MemoryUltraLow,
        AFFINE_P4,
        None
    ),
    reference_case!(
        "test.biwfa.affine.p5",
        AlignmentScope::Alignment,
        MemoryModel::MemoryUltraLow,
        AFFINE_P5,
        None
    ),
    reference_case!(
        "test.biwfa.affine.wfapt0",
        AlignmentScope::Alignment,
        MemoryModel::MemoryUltraLow,
        AFFINE_DEFAULT,
        WFAPT0
    ),
    reference_case!(
        "test.biwfa.affine.wfapt1",
        AlignmentScope::Alignment,
        MemoryModel::MemoryUltraLow,
        AFFINE_DEFAULT,
        WFAPT1
    ),
    reference_case!(
        "test.biwfa.score.indel",
        AlignmentScope::Score,
        MemoryModel::MemoryUltraLow,
        Penalties::Indel,
        None
    ),
    reference_case!(
        "test.biwfa.score.edit",
        AlignmentScope::Score,
        MemoryModel::MemoryUltraLow,
        Penalties::Edit,
        None
    ),
    reference_case!(
        "test.biwfa.score.affine",
        AlignmentScope::Score,
        MemoryModel::MemoryUltraLow,
        AFFINE_DEFAULT,
        None
    ),
    reference_case!(
        "test.biwfa.score.affine2p",
        AlignmentScope::Score,
        MemoryModel::MemoryUltraLow,
        AFFINE2P_DEFAULT,
        None
    ),
    reference_case!(
        "test.biwfa.score.affine.p0",
        AlignmentScope::Score,
        MemoryModel::MemoryUltraLow,
        AFFINE_P0,
        None
    ),
    reference_case!(
        "test.biwfa.score.affine.p1",
        AlignmentScope::Score,
        MemoryModel::MemoryUltraLow,
        AFFINE_P1,
        None
    ),
    reference_case!(
        "test.biwfa.score.affine.p2",
        AlignmentScope::Score,
        MemoryModel::MemoryUltraLow,
        AFFINE_P2,
        None
    ),
    reference_case!(
        "test.biwfa.score.affine.p3",
        AlignmentScope::Score,
        MemoryModel::MemoryUltraLow,
        AFFINE_P3,
        None
    ),
    reference_case!(
        "test.biwfa.score.affine.p4",
        AlignmentScope::Score,
        MemoryModel::MemoryUltraLow,
        AFFINE_P4,
        None
    ),
    reference_case!(
        "test.biwfa.score.affine.p5",
        AlignmentScope::Score,
        MemoryModel::MemoryUltraLow,
        AFFINE_P5,
        None
    ),
    reference_case!(
        "test.biwfa.score.affine.wfapt0",
        AlignmentScope::Score,
        MemoryModel::MemoryUltraLow,
        AFFINE_DEFAULT,
        WFAPT0
    ),
    reference_case!(
        "test.biwfa.score.affine.wfapt1",
        AlignmentScope::Score,
        MemoryModel::MemoryUltraLow,
        AFFINE_DEFAULT,
        WFAPT1
    ),
];

#[test]
#[ignore = "runs the full upstream WFA2 reference matrix"]
fn wfa2_reference_golden_suite() {
    let sequence_pairs = parse_sequence_pairs(SEQUENCES);

    for case in CASES {
        let expected_rows = parse_expected_rows(case.name, case.expected);
        assert_eq!(
            expected_rows.len(),
            sequence_pairs.len(),
            "case {} fixture row count does not match sequence-pair count",
            case.name
        );

        let mut aligner = build_aligner(*case);
        for (row_index, ((pattern, text), expected)) in
            sequence_pairs.iter().zip(expected_rows.iter()).enumerate()
        {
            let result = aligner.align_end_to_end(pattern.as_bytes(), text.as_bytes());
            if result.status != AlignmentStatus::StatusAlgCompleted {
                panic!(
                    "WFA2 reference status mismatch\ncase: {}\nrow: {}\nstatus: {}\npattern: {}\ntext: {}",
                    case.name,
                    row_index + 1,
                    result.status,
                    pattern,
                    text
                );
            }

            let actual_score = aligner.score();
            let actual_cigar = if case.scope == AlignmentScope::Score {
                "-".to_string()
            } else {
                format_wfa_cigar(&aligner.wfa_cigar_bytes())
            };

            if actual_score != expected.score || actual_cigar != expected.cigar {
                panic!(
                    "WFA2 reference mismatch\ncase: {}\nrow: {}\nexpected: {}\t{}\nactual: {}\t{}\npattern: {}\ntext: {}",
                    case.name,
                    row_index + 1,
                    expected.score,
                    expected.cigar,
                    actual_score,
                    actual_cigar,
                    pattern,
                    text
                );
            }
        }
    }
}

fn build_aligner(case: Case) -> WFAligner {
    let builder = WFAligner::builder(case.scope, case.memory);
    let builder = match case.penalties {
        Penalties::Indel => builder.indel(),
        Penalties::Edit => builder.edit(),
        Penalties::Affine {
            match_,
            mismatch,
            gap_opening,
            gap_extension,
        } => builder.affine_with_match(match_, mismatch, gap_opening, gap_extension),
        Penalties::Affine2p {
            match_,
            mismatch,
            gap_opening1,
            gap_extension1,
            gap_opening2,
            gap_extension2,
        } => builder.affine2p_with_match(
            match_,
            mismatch,
            gap_opening1,
            gap_extension1,
            gap_opening2,
            gap_extension2,
        ),
    };

    let builder = match case.heuristic {
        Some(ReferenceHeuristic::WfAdaptive {
            min_wavefront_length,
            max_distance_threshold,
            steps_between_cutoffs,
        }) => builder.with_heuristics(Heuristics::wf_adaptive(
            steps_between_cutoffs,
            min_wavefront_length,
            max_distance_threshold,
        )),
        None => builder,
    };

    builder.build().unwrap_or_else(|err| {
        panic!(
            "failed to build WFA2 reference aligner for case {}: {}",
            case.name, err
        )
    })
}

fn parse_sequence_pairs(input: &'static str) -> Vec<(&'static str, &'static str)> {
    let mut pairs = Vec::new();
    let mut lines = input.lines().enumerate();

    while let Some((pattern_line_number, pattern_line)) = lines.next() {
        let Some(pattern) = pattern_line.strip_prefix('>') else {
            panic!(
                "expected pattern line {} to start with '>'",
                pattern_line_number + 1
            );
        };

        let Some((text_line_number, text_line)) = lines.next() else {
            panic!(
                "missing text line after pattern line {}",
                pattern_line_number + 1
            );
        };
        let Some(text) = text_line.strip_prefix('<') else {
            panic!(
                "expected text line {} to start with '<'",
                text_line_number + 1
            );
        };

        pairs.push((pattern, text));
    }

    pairs
}

fn parse_expected_rows<'a>(case_name: &str, input: &'a str) -> Vec<ExpectedRow<'a>> {
    input
        .lines()
        .enumerate()
        .map(|(index, line)| {
            let (score, cigar) = line.split_once('\t').unwrap_or_else(|| {
                panic!(
                    "expected tab-separated score and CIGAR in case {} row {}",
                    case_name,
                    index + 1
                )
            });
            let score = score.parse::<i32>().unwrap_or_else(|err| {
                panic!(
                    "invalid score in case {} row {}: {}",
                    case_name,
                    index + 1,
                    err
                )
            });
            ExpectedRow { score, cigar }
        })
        .collect()
}

fn format_wfa_cigar(operations: &[u8]) -> String {
    let Some((&first, rest)) = operations.split_first() else {
        return "-".to_string();
    };

    let mut cigar = String::new();
    let mut current = first;
    let mut length = 1usize;

    for &operation in rest {
        if operation == current {
            length += 1;
        } else {
            push_cigar_run(&mut cigar, current, length);
            current = operation;
            length = 1;
        }
    }

    push_cigar_run(&mut cigar, current, length);
    cigar
}

fn push_cigar_run(cigar: &mut String, operation: u8, length: usize) {
    match operation {
        b'M' | b'X' | b'I' | b'D' => {
            cigar.push_str(&length.to_string());
            cigar.push(operation as char);
        }
        _ => panic!("unexpected WFA CIGAR operation byte: {}", operation),
    }
}
