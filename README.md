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
```
