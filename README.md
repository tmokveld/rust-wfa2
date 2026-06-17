# Rust bindings for WFA2-Lib

Rust language bindings for the excellent
[WFA2-Lib](https://github.com/smarco/WFA2-lib) library.

Work in progress. Tests and features are not yet complete.

## Autovectorization

Remember to specify the correct C compiler! For me it is `CC=/usr/local/opt/llvm/bin/clang`.

## OpenMP

OpenMP support is opt-in. The default build stays serial and avoids an OpenMP
runtime dependency.

```sh
cargo build --release --features openmp
```

In testing, WFA2's OpenMP path did not provide reliable speedups for the
workloads tried.

On Linux, GCC/libgomp is the default OpenMP runtime:

```sh
CC=gcc CXX=g++ cargo build --release --features openmp
```

On macOS with Homebrew LLVM/libomp:

```sh
LLVM_PREFIX="$(brew --prefix llvm)" \
LIBOMP_PREFIX="$(brew --prefix libomp)" \
CC="$(brew --prefix llvm)/bin/clang" \
CXX="$(brew --prefix llvm)/bin/clang++" \
cargo build --release --features openmp
```

Use `WFA2_OPENMP_LIB=omp` or `WFA2_OPENMP_LIB=gomp` to override the runtime
linked by Cargo, and `WFA2_OPENMP_LIB_DIR` if the runtime is in a non-standard
library directory.

## Usage

```rust
use rust_wfa2::aligner::{AlignmentScope, AlignmentStatus, MemoryModel, WFAligner};

let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
    .affine(6, 4, 2)
    .build().unwrap();

let pattern = b"TCTTTACTCGCGCGTTGGAGAAATACAATAGT";
let text = b"TCTATACTGCGCGTTTGGAGAAATAAAATAGT";
let result = aligner.align_end_to_end(pattern, text);
assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
assert_eq!(aligner.score(), -24);
assert_eq!(
    aligner.cigar_operations(),
    b"MMMXMMMMDMMMMMMMIMMMMMMMMMXMMMMMM"
);
```

### CIGAR orientation and SAM

WFA2 CIGAR operations describe how to transform the `pattern` argument into the
`text` argument. In this Rust wrapper, `pattern` is usually the query and `text`
is usually the reference, so raw WFA2 operations use the opposite insertion and
deletion orientation from SAM's reference-to-query CIGAR semantics.

`get_sam_cigar()` uses BAM/SAM's packed integer encoding, but it does not change
that WFA2 operation orientation. For SAM-compliant reference/query CIGAR output,
either call the aligner with `pattern = reference` and `text = query`, or swap
`I` and `D` after decoding.

### WFA2 plot dumps

WFA2's native wavefront plot recorder is available for debugging and tooling.
It writes WFA2's `.plot` text format; PNG rendering is left to external tools
such as WFA2's `scripts/wfa.plot.py`.

```rust
use rust_wfa2::aligner::{AlignmentScope, AlignmentStatus, MemoryModel, PlotOptions, WFAligner};

let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
    .with_plotting(PlotOptions::default())
    .affine(6, 4, 2)
    .build().unwrap();

let result = aligner.align_end_to_end(
    b"TCTTTACTCGCGCGTTGGAGAAATACAATAGT",
    b"TCTATACTGCGCGTTTGGAGAAATAAAATAGT",
);
assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
aligner.write_plot("debug.plot").unwrap();
```

Heuristics are opt-in. WFA2 supports combining at most one adaptive heuristic,
one drop heuristic, and one band heuristic:

```rust
use rust_wfa2::aligner::{
    AdaptiveHeuristic, AlignmentScope, BandHeuristic, DropHeuristic, Heuristics, MemoryModel,
    WFAligner,
};

let heuristics = Heuristics::new(10)
    .with_adaptive(AdaptiveHeuristic::WfAdaptive {
        min_wavefront_length: 10,
        max_distance_threshold: 50,
    })
    .with_drop(DropHeuristic::XDrop { xdrop: 100 })
    .with_band(BandHeuristic::Static {
        min_k: -50,
        max_k: 50,
    });

let mut aligner = WFAligner::builder(AlignmentScope::Alignment, MemoryModel::MemoryLow)
    .affine(6, 4, 2)
    .with_heuristics(heuristics)
    .build().unwrap();
```
