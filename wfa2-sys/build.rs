use std::collections::HashSet;
use std::env;
use std::path::{Path, PathBuf};

#[derive(Debug)]
struct IgnoreMacros(HashSet<String>);

impl bindgen::callbacks::ParseCallbacks for IgnoreMacros {
    fn will_parse_macro(&self, name: &str) -> bindgen::callbacks::MacroParsingBehavior {
        if self.0.contains(name) {
            bindgen::callbacks::MacroParsingBehavior::Ignore
        } else {
            bindgen::callbacks::MacroParsingBehavior::Default
        }
    }
}

fn main() {
    let openmp_enabled = env::var_os("CARGO_FEATURE_OPENMP").is_some();
    emit_rerun_if_env_changed();
    println!("cargo:rerun-if-changed=WFA2-lib");

    let mut cmake = cmake::Config::new("WFA2-lib");
    // Force an optimized C build regardless of the Cargo profile. WFA2's CMakeLists only appends
    // OPTIMIZE_FLAGS (including EXTRA_FLAGS) for Release/RelWithDebInfo builds.
    cmake.profile("Release");

    // EXTRA_FLAGS is a cmake cache variable consumed by WFA2's CMakeLists, not a C preprocessor
    // flag. It must be passed via `define` so cmake appends it to CMAKE_C_FLAGS/CMAKE_CXX_FLAGS;
    // passing it through `cflag` only defines an unused `-DEXTRA_FLAGS=...` macro.
    let mut extra_flags = String::from("-ftree-vectorize");
    if target_arch() == "x86_64" {
        // SSE tuning is x86-only. Enabling it on aarch64 (e.g. Apple Silicon) fails the build.
        extra_flags.push_str(" -msse2 -mfpmath=sse");
    }
    for flag in native_codegen_flags() {
        extra_flags.push(' ');
        extra_flags.push_str(flag);
    }
    cmake.define("EXTRA_FLAGS", &extra_flags);

    if openmp_enabled {
        cmake.define("OPENMP", "ON");
        configure_openmp_cmake(&mut cmake);
    }

    let out_dir = cmake.build();
    println!("cargo:rustc-link-search=native={}/lib", out_dir.display());
    println!("cargo:rustc-link-search=native={}/lib64", out_dir.display());
    println!("cargo:rustc-link-lib=static=wfa2");
    if openmp_enabled {
        emit_openmp_linking();
    }

    let ignored_macros = IgnoreMacros(
        vec![
            "FP_INFINITE".into(),
            "FP_NAN".into(),
            "FP_NORMAL".into(),
            "FP_SUBNORMAL".into(),
            "FP_ZERO".into(),
            "IPPORT_RESERVED".into(),
        ]
        .into_iter()
        .collect(),
    );

    bindgen::Builder::default()
        .header("WFA2-lib/utils/commons.h")
        .header("WFA2-lib/wavefront/wfa.h")
        .clang_arg("--include-directory=WFA2-lib")
        .parse_callbacks(Box::new(ignored_macros))
        .generate()
        .expect("Unable to generate bindings")
        .write_to_file(out_dir.join("bindings.rs"))
        .expect("Couldn't write bindings!");
}

fn emit_rerun_if_env_changed() {
    for var in [
        "CC",
        "CXX",
        "CMAKE_PREFIX_PATH",
        "HOMEBREW_PREFIX",
        "LLVM_PREFIX",
        "LIBOMP_PREFIX",
        "WFA2_OPENMP_LIB",
        "WFA2_OPENMP_LIB_DIR",
    ] {
        println!("cargo:rerun-if-env-changed={var}");
    }
}

fn configure_openmp_cmake(cmake: &mut cmake::Config) {
    if target_os() != "macos" {
        return;
    }

    if env::var_os("CC").is_none() {
        if let Some(clang) = prefixed_file("LLVM_PREFIX", "llvm", "bin/clang") {
            cmake.define("CMAKE_C_COMPILER", clang);
        }
    }

    if env::var_os("CXX").is_none() {
        if let Some(clangxx) = prefixed_file("LLVM_PREFIX", "llvm", "bin/clang++") {
            cmake.define("CMAKE_CXX_COMPILER", clangxx);
        }
    }

    if env::var_os("CMAKE_PREFIX_PATH").is_none() {
        let mut prefix_paths = Vec::new();
        if let Some(llvm_prefix) = package_prefix("LLVM_PREFIX", "llvm") {
            prefix_paths.push(llvm_prefix);
        }
        if let Some(libomp_prefix) = package_prefix("LIBOMP_PREFIX", "libomp") {
            prefix_paths.push(libomp_prefix);
        }

        if !prefix_paths.is_empty() {
            let prefix_path = prefix_paths
                .iter()
                .map(|path| path.to_string_lossy())
                .collect::<Vec<_>>()
                .join(";");
            cmake.define("CMAKE_PREFIX_PATH", prefix_path);
        }
    }
}

fn emit_openmp_linking() {
    if let Some(lib_dir) = env::var_os("WFA2_OPENMP_LIB_DIR") {
        println!(
            "cargo:rustc-link-search=native={}",
            PathBuf::from(lib_dir).display()
        );
    }

    let target_os = target_os();
    if target_os == "macos" {
        if let Some(libomp_prefix) = package_prefix("LIBOMP_PREFIX", "libomp") {
            println!(
                "cargo:rustc-link-search=native={}",
                libomp_prefix.join("lib").display()
            );
        }
    }

    let openmp_lib = env::var("WFA2_OPENMP_LIB").unwrap_or_else(|_| default_openmp_lib());
    println!("cargo:rustc-link-lib=dylib={openmp_lib}");
}

fn default_openmp_lib() -> String {
    match target_os().as_str() {
        "macos" => "omp".to_string(),
        "linux" if compiler_looks_like_clang() => "omp".to_string(),
        "linux" => "gomp".to_string(),
        _ => "omp".to_string(),
    }
}

fn compiler_looks_like_clang() -> bool {
    env::var("CC")
        .map(|cc| cc.to_ascii_lowercase().contains("clang"))
        .unwrap_or(false)
}

fn target_os() -> String {
    env::var("CARGO_CFG_TARGET_OS").unwrap_or_default()
}

fn target_arch() -> String {
    env::var("CARGO_CFG_TARGET_ARCH").unwrap_or_default()
}

/// Opt-in native codegen flags for WFA2's SIMD extend kernels (off by default).
///
/// `native` (`-march=native`) tunes for the building machine and is the superset, so it
/// takes precedence over `avx2`. `avx2` (`-mavx2`) only makes sense on x86_64, where it
/// defines `__AVX2__` and activates WFA2's hand-written AVX2 kernel; on other targets it is
/// ignored with a warning rather than failing the build.
fn native_codegen_flags() -> Vec<&'static str> {
    let native = env::var_os("CARGO_FEATURE_NATIVE").is_some();
    let avx2 = env::var_os("CARGO_FEATURE_AVX2").is_some();

    if native {
        return vec!["-march=native"];
    }
    if avx2 {
        if target_arch() == "x86_64" {
            return vec!["-mavx2"];
        }
        println!(
            "cargo:warning=wfa2-sys: the `avx2` feature only affects x86_64 targets; ignoring on {}",
            target_arch()
        );
    }
    Vec::new()
}

fn prefixed_file(env_var: &str, package: &str, relative_path: &str) -> Option<PathBuf> {
    let path = package_prefix(env_var, package)?.join(relative_path);
    path.is_file().then_some(path)
}

fn package_prefix(env_var: &str, package: &str) -> Option<PathBuf> {
    if let Some(prefix) = env_path(env_var) {
        return Some(prefix);
    }

    homebrew_package_prefix(package)
}

fn env_path(env_var: &str) -> Option<PathBuf> {
    env::var_os(env_var)
        .map(PathBuf::from)
        .filter(|path| path.exists())
}

fn homebrew_package_prefix(package: &str) -> Option<PathBuf> {
    let mut candidates = Vec::new();
    if let Some(homebrew_prefix) = env_path("HOMEBREW_PREFIX") {
        candidates.push(homebrew_prefix.join("opt").join(package));
    }
    candidates.push(Path::new("/opt/homebrew/opt").join(package));
    candidates.push(Path::new("/usr/local/opt").join(package));

    candidates.into_iter().find(|path| path.exists())
}
