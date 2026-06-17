# Rust bindings for WFA2-Lib

Rust language bindings for the excellent
[WFA2-Lib](https://github.com/smarco/WFA2-lib) library.

Work in progress. Tests and features are not yet complete.

## Native codegen (SIMD)

The WFA2 C library is always compiled in an optimized (`-O3`,
release) configuration with autovectorization enabled, so the default build is
portable and runs on any CPU of the target architecture. No extra flags are
required.

WFA2 also ships a AVX2 extend kernel that is only compiled in when
the compiler defines `__AVX2__` (i.e. when AVX2 codegen is enabled). Because
emitting AVX2 instructions makes the binary crash on CPUs without AVX2, this is
**opt-in** via two features (both off by default):

```sh
# Portable to any AVX2-capable x86_64 CPU (Haswell, 2013+):
cargo build --release --features avx2

# Tune for the building machine specifically (-march=native).
# Fastest, but the resulting binary is not portable to other CPUs:
cargo build --release --features native
```

Notes:

- `avx2` only affects `x86_64`; on other architectures it is ignored with a
  build warning.
- `native` works on any architecture and is the superset, so it takes
  precedence if both are enabled.
- You can still override the C compiler if needed, e.g.
  `CC=/usr/local/opt/llvm/bin/clang`.

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

let query = b"TCTTTACTCGCGCGTTGGAGAAATACAATAGT";
let reference = b"TCTATACTGCGCGTTTGGAGAAATAAAATAGT";
let result = aligner.align_end_to_end(query, reference);
assert_eq!(result.status, AlignmentStatus::StatusAlgCompleted);
assert_eq!(aligner.score(), -24);
assert_eq!(
    aligner.sam_cigar_bytes(),
    b"MMMXMMMMIMMMMMMMMDMMMMMMMMMXMMMMMM"
);
```

### CIGAR orientation and SAM

WFA2 CIGAR operations describe how to transform the `pattern` argument into the
`text` argument. In this Rust wrapper, `pattern` is usually the query and `text`
is usually the reference, so raw WFA2 operations use the opposite insertion and
deletion orientation from SAM's reference-to-query CIGAR semantics.

Use `sam_cigar()` or `sam_packed_cigar()` when the last alignment was called with `pattern = query` and `text = reference`. Use `wfa_cigar_bytes()`,
`wfa_cigar()`, or `wfa_packed_cigar()` when you want WFA's native
pattern-to-text orientation, including reference-first workflows.

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

Heuristics are opt-in. This intentionally differs from WFA2's C
`wavefront_aligner_attr_default`, which enables WF-adaptive pruning with
`min_wavefront_length = 10`, `max_distance_threshold = 50`, and
`steps_between_cutoffs = 1`. The Rust wrapper defaults to exact alignment unless
you configure heuristics explicitly; use `Heuristics::wfa2_default()` to recover
the C default behavior.

WFA2 supports combining at most one adaptive heuristic, one drop heuristic, and
one band heuristic:

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
