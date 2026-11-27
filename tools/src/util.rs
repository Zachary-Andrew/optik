//! Small process-level helpers.
//! They choose default thread counts, set the OpenMP budget, track peak memory, parse sizes such as 8G and detect when two paths name the same file.

use super::*;

/// Makes the OpenMP runtime follow the thread budget `threads`.
///
/// libgomp reads OMP_NUM_THREADS and OMP_WAIT_POLICY once, in its load-time constructor before main, so setting them here would change nothing.
/// When they do not already say what -t needs, the process re-executes itself once with them set, keeping the same arguments, standard streams and the rest of the environment.
/// OMP_NUM_THREADS = t becomes the default of a region entered from any thread.
/// OMP_WAIT_POLICY = PASSIVE (unless the user set one) makes idle OpenMP threads sleep immediately instead of spinning.
/// By default they spin after a region, and a phase starting its own threads then ran beside them, past the budget, as bench/thread_check.py found.
/// OPTIK_OMP_ENV marks the re-executed process, so it never loops.
pub(super) fn omp_env_for_budget(threads: usize) {
    use std::os::unix::process::CommandExt;
    let want_n = threads.to_string();
    let n_ok = std::env::var("OMP_NUM_THREADS").ok().as_deref() == Some(want_n.as_str());
    let wait_ok = std::env::var_os("OMP_WAIT_POLICY").is_some();

    if (n_ok && wait_ok) || std::env::var_os("OPTIK_OMP_ENV").is_some() {
        return;
    }

    let Ok(exe) = std::env::current_exe() else {
        return;
    };
    let mut cmd = std::process::Command::new(exe);

    cmd.args(std::env::args_os().skip(1))
        .env("OMP_NUM_THREADS", &want_n)
        .env("OPTIK_OMP_ENV", "1");

    if !wait_ok {
        cmd.env("OMP_WAIT_POLICY", "PASSIVE");
    }

    // exec replaces this process under the same pid and returns only on failure.
    // The run then goes on as it is, with omp_set_num_threads still covering this thread.
    let err = cmd.exec();

    eprintln!(
        "[optik] warning: could not re-execute with the OpenMP environment ({err}); continuing"
    );
}

/// Default thread budget, half the logical CPUs and at least 1.
/// -t is a budget for everything a command keeps running at once, i.e. worker pools, input readers (a .gz is decoded inside its reader's thread, with no helper process) and the C++ core's OpenMP regions, which bench/thread_check.py measures.
/// The C++ core resolves 0 the same way (thread_budget.hpp).
pub(super) fn default_threads() -> usize {
    // One thread per two logical CPUs, i.e. per physical core on an SMT machine.
    const LOGICAL_CPUS_PER_THREAD: usize = 2;
    let cpus = std::thread::available_parallelism()
        .map(|count| count.get())
        .unwrap_or(1);

    (cpus / LOGICAL_CPUS_PER_THREAD).max(1)
}

/// Peak resident-set growth since creation.
///
/// Reads VmRSS from /proc/self/status at creation and at each `sample()`, so the difference is the memory an operation added, without the process's baseline (libraries, runtime).
pub(crate) struct RssWatermark {
    pub(super) before_kb: u64,
    pub(super) peak_kb: u64,
}

impl RssWatermark {
    /// Starts tracking from the current resident set.
    pub(crate) fn new() -> Self {
        let rss_kb = Self::read_kb();

        RssWatermark {
            before_kb: rss_kb,
            peak_kb: rss_kb,
        }
    }

    /// Records the current RSS, updating the peak if it increased.
    pub(crate) fn sample(&mut self) {
        let rss_kb = Self::read_kb();

        if rss_kb > self.peak_kb {
            self.peak_kb = rss_kb;
        }
    }

    /// Net memory increase from start to peak, in megabytes.
    pub(crate) fn net_mb(&self) -> f64 {
        self.peak_kb.saturating_sub(self.before_kb) as f64 / KB_PER_MB
    }

    /// Reads the current RSS from the kernel's proc filesystem.
    pub(super) fn read_kb() -> u64 {
        std::fs::read_to_string("/proc/self/status")
            .unwrap_or_default()
            .lines()
            .find(|line| line.starts_with("VmRSS:"))
            .and_then(|line| line.split_whitespace().nth(1))
            .and_then(|kb| kb.parse().ok())
            .unwrap_or(0)
    }
}

/// Parses a memory size, a number with an optional K/M/G/T suffix (powers of 1024, where a bare number is bytes).
/// Returns `None` if the size is malformed or not positive.
pub(super) fn parse_size(input: &str) -> Option<u64> {
    // 8G, 8GB, 8GiB and 8g are the same size.
    let trimmed = input.trim();
    let without_ib = trimmed
        .strip_suffix("iB")
        .or_else(|| trimmed.strip_suffix("ib"))
        .unwrap_or(trimmed);
    // A trailing B counts as a unit only after a prefix letter. A bare 512B stays invalid.
    let number_unit = if without_ib.len() > 1
        && without_ib[..without_ib.len() - 1].ends_with(|ch: char| "kmgtKMGT".contains(ch))
    {
        without_ib
            .strip_suffix('B')
            .or_else(|| without_ib.strip_suffix('b'))
            .unwrap_or(without_ib)
    } else {
        without_ib
    };
    let (num, mul) = match number_unit.chars().last()?.to_ascii_uppercase() {
        'K' => (&number_unit[..number_unit.len() - 1], KIB),
        'M' => (&number_unit[..number_unit.len() - 1], MIB),
        'G' => (&number_unit[..number_unit.len() - 1], GIB),
        'T' => (&number_unit[..number_unit.len() - 1], TIB),
        _ => (number_unit, 1u64),
    };
    let value: f64 = num.parse().ok()?;

    if !value.is_finite() || value <= 0.0 {
        return None;
    }

    let bytes = value * mul as f64;

    // The f64 cast saturates silently, so overflow is caught here.
    if bytes >= u64::MAX as f64 {
        return None;
    }

    Some(bytes as u64)
}

/// Whether two paths name the same existing file, by device and inode.
pub(super) fn same_file(path_a: &str, path_b: &str) -> bool {
    use std::os::unix::fs::MetadataExt;

    match (std::fs::metadata(path_a), std::fs::metadata(path_b)) {
        (Ok(meta_a), Ok(meta_b)) => meta_a.dev() == meta_b.dev() && meta_a.ino() == meta_b.ino(),
        _ => false,
    }
}
