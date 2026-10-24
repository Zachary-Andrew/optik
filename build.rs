//! Build script for optik's native C++ code.
//!
//! Compiles src/*.cpp with the cc crate, reads the SIMD capability of the build machine from /proc/cpuinfo, on Linux only, with the scalar path used elsewhere.
//! It enables AVX-512 F+BW+DQ or AVX2 accordingly, then links libgomp and zlib, the zlib-ng in temp/zlib-ng when it is built.
fn main() {
    // Nix's compiler wrapper drops -march=native unless NIX_ENFORCE_NO_NATIVE is 0, which would build the kernels without this CPU's instruction set.
    // The variable has no effect outside Nix.
    // SAFETY: the build script is single-threaded here, before cc starts any compiler.
    unsafe { std::env::set_var("NIX_ENFORCE_NO_NATIVE", "0") };
    println!("cargo:rerun-if-env-changed=NIX_ENFORCE_NO_NATIVE");
    // Re-run if any header changes, since cc caches by .cpp mtime only
    for dir in ["include", "src"] {
        if let Ok(rd) = std::fs::read_dir(dir) {
            for e in rd.flatten() {
                let p = e.path();
                if p.extension().is_some_and(|x| x == "hpp" || x == "cpp") {
                    println!("cargo:rerun-if-changed={}", p.display());
                }
            }
        }
    }
    // Read actual CPU flags from /proc/cpuinfo once and derive SIMD tiers.
    // The C++ code guards every intrinsic behind its preprocessor macro, so the scalar fallback is always compiled in and the binary is correct on any x86-64 CPU.
    let flags = cpu_flags();

    let avx512f = flags.contains("avx512f");
    let avx512bw = flags.contains("avx512bw");
    let avx512dq = flags.contains("avx512dq");
    let avx2 = flags.contains("avx2");

    let zlib = zlib_location();
    probe_dependencies(zlib.as_ref().map(|(dir, _)| dir.as_str()));

    let mut build = cc::Build::new();
    build
        .cpp(true)
        .file("src/ffi.cpp")
        .file("src/kmer_count.cpp")
        .file("src/ffi_dbg.cpp")
        .file("src/ffi_fm.cpp")
        .file("src/ffi_map.cpp")
        .file("src/ffi_spss.cpp")
        .include("include")
        .flag("-std=c++17")
        .flag("-O3")
        .flag("-march=native")
        .flag("-funroll-loops")
        .flag("-Wall");
    if let Some((dir, _)) = &zlib {
        build.include(format!("{dir}/include"));
    }

    // OpenMP for the parallel build phases.
    build.flag("-fopenmp");

    // Enable SIMD tiers only when the host CPU supports them.
    // The C++ headers check __AVX512F__ / __AVX2__ before using any intrinsic, so adding a flag here is both necessary and sufficient.
    if avx512f && avx512bw && avx512dq {
        build
            .flag("-mavx512f")
            .flag("-mavx512bw")
            .flag("-mavx512dq");
        println!("cargo:warning=optik: AVX-512 (F+BW+DQ) enabled");
    } else if avx512f && avx512bw {
        build.flag("-mavx512f").flag("-mavx512bw");
        println!("cargo:warning=optik: AVX-512 (F+BW) enabled, no DQ");
    } else if avx2 {
        build.flag("-mavx2");
        println!("cargo:warning=optik: AVX2 enabled");
    } else {
        println!("cargo:warning=optik: no SIMD extensions, using scalar path");
    }

    build.compile("optik_native");

    // Link libgomp after liboptik_native.
    // rustc release links pass --as-needed, which would otherwise drop libgomp (no omp references are visible before the archive that introduces them is scanned) and leave ffi.o's omp_* symbols unresolved.
    println!("cargo:rustc-link-lib=gomp");
    // zlib decodes .gz input in the reader's own thread (fasta_reader.hpp), so a reader is one thread of the -t budget.
    // zlib_location() picks the library, and the system libz is linked when it finds none.
    match &zlib {
        Some((dir, lib)) => {
            println!("cargo:rustc-link-search=native={lib}");
            if dir.ends_with("temp/zlib-ng") {
                println!("cargo:rustc-link-lib=static=z");
                println!("cargo:warning=optik: zlib-ng from {dir}");
            } else {
                println!("cargo:rustc-link-lib=z");
                println!("cargo:warning=optik: zlib from {dir}");
            }
        }
        None => println!("cargo:rustc-link-lib=z"),
    }
}

/// The headers optik's C++ code needs, each with the packages that provide it
const PROBES: [(&str, &str, &str); 3] = [
    (
        "c++",
        "cstdint",
        "a C++17 compiler: build-essential (apt), gcc-c++ (dnf, zypper), base-devel (pacman)",
    ),
    (
        "openmp",
        "omp.h",
        "OpenMP: GCC's libgomp, part of gcc on most distributions (dnf: libgomp; clang needs libomp-dev)",
    ),
    (
        "zlib",
        "zlib.h",
        "zlib headers: zlib1g-dev (apt), zlib-devel (dnf, zypper), zlib (pacman), or set ZLIB_DIR",
    ),
];

/// probe_dependencies() - Stop the build early when a compiler, OpenMP or zlib is missing.
/// @zlib_dir: the zlib prefix from zlib_location(), or None for the system zlib.
///
/// Each probe compiles a one-line file that includes one header, so a missing package is named with the command that installs it, instead of surfacing as a compiler or linker error deep in the build.
/// The README's Install section has the full commands per platform.
fn probe_dependencies(zlib_dir: Option<&str>) {
    let out = std::env::var("OUT_DIR").expect("cargo sets OUT_DIR");
    let mut missing = Vec::new();
    for (name, header, hint) in PROBES {
        let src = format!("{out}/probe_{name}.cpp");
        let body = format!("#include <{header}>\nint optik_probe() {{ return 0; }}\n");
        if std::fs::write(&src, body).is_err() {
            return;
        }
        let mut b = cc::Build::new();
        b.cpp(true)
            .flag("-std=c++17")
            .flag("-fopenmp")
            .cargo_metadata(false)
            .warnings(false)
            .file(&src);
        if let Some(dir) = zlib_dir {
            b.include(format!("{dir}/include"));
        }
        if b.try_compile(&format!("optik_probe_{}", name.replace('+', "x")))
            .is_err()
        {
            missing.push(hint);
            if name == "c++" {
                break; // without a compiler the other probes say nothing new
            }
        }
    }
    if missing.is_empty() {
        return;
    }
    for hint in &missing {
        println!("cargo:warning=optik: missing {hint}");
    }
    eprintln!(
        "\noptik cannot be built until these are installed:\n{}\n\nSee README.md, section Install, \
         for the command on your platform.\n",
        missing
            .iter()
            .map(|m| format!("  - {m}"))
            .collect::<Vec<_>>()
            .join("\n")
    );
    std::process::exit(1);
}

/// The zlib to build against, as (install directory, library directory).
///
/// ZLIB_DIR (include/ and lib/) wins when set, and any zlib-compatible library works there.
/// Otherwise the zlib-ng that `bench/build_tools.sh zlib-ng` installs in temp/zlib-ng is used when present, a faster inflate in zlib-compatible mode.
/// None means the system zlib.
fn zlib_location() -> Option<(String, String)> {
    println!("cargo:rerun-if-env-changed=ZLIB_DIR");
    if let Ok(d) = std::env::var("ZLIB_DIR") {
        return Some((d.clone(), format!("{d}/lib")));
    }
    let root = std::env::var("CARGO_MANIFEST_DIR").ok()?;
    let dir = format!("{root}/temp/zlib-ng");
    println!("cargo:rerun-if-changed={dir}");
    ["lib", "lib64"]
        .iter()
        .map(|l| format!("{dir}/{l}"))
        .find(|l| std::path::Path::new(&format!("{l}/libz.a")).exists())
        .map(|lib| (dir.clone(), lib))
}

/// The contents of /proc/cpuinfo, searched for flag names, or an empty string where it does not exist.
fn cpu_flags() -> String {
    // /proc/cpuinfo is Linux-only.  On other platforms return empty string
    // so every SIMD tier stays disabled and the scalar path is used.
    std::fs::read_to_string("/proc/cpuinfo").unwrap_or_default()
}
