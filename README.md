# Optik

Optik is a genome analysis toolkit for k-mer indexing, de Bruijn graph construction and read mapping, built on cache line open-addressing hash tables and a minimal perfect hash whose per-bucket pilots are inspired by tiny pointers.

| Command | What it does |
|---|---|
| `build-spss`, `query` | Exact k-mer membership and dense k-mer ids |
| `build-fm`, `query` | The same membership answers from a smaller but slower FM index |
| `spss-update` | Adds or removes k-mers in an SPSS index without the original input |
| `dbruijn` | Counts k-mers and writes the compacted de Bruijn graph as GFA or unitig FASTA |
| `map` | Aligns short, paired-end and long reads, writing SAM |
| `tui` | A terminal interface over all of the above |

- [Install](#install)
- [Quick start](#quick-start) and [control flags](#control-flags)
- [Usage](#usage) for [k-mer indexes](#k-mer-indexes), [de Bruijn graphs](#de-bruijn-graphs) and [read mapping](#read-mapping)
- [Intended use and limitations](#intended-use-and-limitations)
- [Benchmarks and validation](#benchmarks-and-validation)

## Install

Optik is compiled from source on Linux x86-64 and needs Rust 1.88 or newer, a C++17 compiler with OpenMP, CMake and zlib.
The packages below provide them.

| Platform | Packages |
|---|---|
| Debian, Ubuntu | `sudo apt install build-essential cmake zlib1g-dev git curl` |
| Fedora, RHEL 10 | `sudo dnf install gcc-c++ make cmake zlib-ng-compat-devel libgomp git curl` |
| RHEL 9 | `sudo dnf install gcc-c++ make cmake zlib-devel libgomp git curl` |
| Arch | `sudo pacman -S --needed base-devel cmake zlib-ng-compat git curl` |
| openSUSE Tumbleweed | `sudo zypper install gcc-c++ make cmake zlib-ng-compat-devel git curl` |
| openSUSE Leap | `sudo zypper install gcc-c++ make cmake zlib-devel git curl` |
| Nix | `nix-shell -p rustc cargo gcc cmake git '(zlib-ng.override { withZlibCompat = true; })'` |

If you do not have rust or a rust package manager the rustup installer below sets up the stable toolchain.

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
. "$HOME/.cargo/env"
```

If rustup came from your package manager instead, `rustup default stable` installs and selects that toolchain.

The Fedora, RHEL 10, Arch, openSUSE Tumbleweed and Nix rows install zlib-ng in zlib-compatible mode, which decodes gzipped input faster than classic zlib with identical results.
With zlib or zlib-ng, build optik with:

```bash
git clone https://github.com/Zachary-Andrew/optik.git && cd optik
cargo build --release          # the binary is target/release/optik
export PATH=$PWD/target/release:$PATH
```

Debian, Ubuntu, RHEL 9 and openSUSE Leap have no zlib-ng package and optik links the system zlib from the installation step, for the speed advantages on those distributions you can use the build script that compiles zlib-ng into `temp/zlib-ng` first that optik then links statically.

```bash
git clone https://github.com/Zachary-Andrew/optik.git && cd optik
bench/build_tools.sh zlib-ng   # zlib-ng into temp/zlib-ng
cargo build --release          # the binary is target/release/optik
export PATH=$PWD/target/release:$PATH
```

If the compiler, OpenMP or zlib is missing, the build stops and names the package to install.

On NixOS and other Nix setups, the compiler wrapper silently drops `-march=native` unless `NIX_ENFORCE_NO_NATIVE=0` is set.
Optik's build and the bench scripts set it for themselves, but other software compiled in the same shell stays generic, so for native builds in general add it to your configuration, for example `environment.variables.NIX_ENFORCE_NO_NATIVE = "0";` in `configuration.nix`.

## Quick start

```bash
optik dbruijn -i genome.fa -k 31 -c 1 -o unitigs.fa         # compacted de Bruijn graph as unitigs
optik build-spss -i unitigs.fa -k 31 -o genome.spss         # index of the genome's 31-mers
optik query -i genome.spss -q reads.fq.gz                   # how many read k-mers occur in the genome
optik map -r genome.fa -q reads.fq.gz -o reads.sam          # the reads aligned to the genome
```

`optik <command> --help` shows each command's flags with examples.
Running `optik` with no arguments in a terminal opens the TUI instead, where `?` or F1 lists the keys.


### Control flags

- Inputs are FASTA or FASTQ, plain or gzipped, and `-` reads standard input.
  Repeating `-i` concatenates files.
- k is 1 to 32, or 1 to 31 for `dbruijn`.
- `-t N` caps how many threads optik keeps busy, and by default it uses half the threads available to it.

## Usage

### k-mer indexes

```bash
optik build-spss -i unitigs.fa -k 31 -o idx.spss                     # from unitigs or any sequences
optik build-spss --compact -i R1.fq.gz -i R2.fq.gz -c 2 -o idx.spss  # count reads, drop k-mers seen once
optik build-fm -i unitigs.fa -k 31 -o idx.fm                         # smaller index, slower lookups
optik query -i idx.spss -q reads.fq                                  # found and missing k-mer counts
optik query -i idx.spss -q reads.fq --ids ids.txt                    # plus a dense id per k-mer window, slower
optik spss-update -i idx.spss --add new.fa --remove old.fa -o idx.spss
```

A k-mer and its reverse complement count as one.
A plain `query` is the fast path, while `--ids` adds an id lookup per window that costs extra time.
`query` takes k from the index and verifies the index's checksum when it loads.
To give each distinct k-mer exactly one id, build with `--compact` or index unitigs (`dbruijn -o x.fa`).

### de Bruijn graphs

```bash
optik dbruijn -i genome.fa -k 31 -c 1 -o graph.gfa        # a genome, keeping every k-mer, as GFA with links
optik dbruijn -i genome.fa -k 31 -c 1 -o unitigs.fa       # the same graph as unitig FASTA
optik dbruijn -i R1.fq.gz -i R2.fq.gz -k 31 -o graph.gfa  # reads, where the default -c 2 drops error k-mers
```

An output named `.fa`, `.fasta` or `.fna` gets unitig FASTA with BCALM-style headers, and any other name gets GFA1.
`--max-mem` limits the k-mer counter, which spills to `--tmp-dir` beyond that limit and otherwise plans with half the available memory.
Neither the thread count nor the memory budget changes the graph.

### Read mapping

```bash
optik map -r ref.fa -q reads.fq -o out.sam                 # short or long reads
optik map -r ref.fa -q R1.fq.gz -q R2.fq.gz -o out.sam     # paired-end (FR)
optik map -r ref.fa -q R1.fq -q R2.fq -R '@RG\tID:s1\tSM:s1' --cs -o out.sam
```

Records carry `NM` and `MD`, `--cs` adds minimap2's `cs` tag, and a split long read gets supplementary records.
Unless `-k` sets it, the seed length is chosen from the reads, 21 for most and 15 for noisy long reads.

## Intended use and limitations

optik is built for k-mer membership over genomes and pangenomes, graphs of reference genomes, and mapping where a confident placement has to be right.
In our benchmarks, which run each compared tool at the settings of its paper, it shows three main strengths.

- SPSS membership queries beat SSHash's on large sets at a similar index size.
- On assembled genomes, `dbruijn` builds the graph 1.7 to 4 times as fast as Cuttlefish.
- At MAPQ ≥ 30, `map` misplaces far fewer simulated human reads than minimap2 or BWA-MEM.

For the following uses, another tool is the better choice.

- **k above 32, or colored, weighted or abundance-annotated sets.**
  optik indexes plain k-mer sets with k up to 32.
- **Dense ids at full speed.**
  `query --ids` returns an id for every window, but those lookups run slower.
- **The smallest possible index.**
  FMSI's index is smaller than `build-fm`'s.
- **Graphs of large sequencing read sets.**
  The genome speedup does not carry over to raw reads, where Cuttlefish matches or beats `dbruijn`'s time on real Illumina sets and needs less memory.
- **Variant calling in segmental duplications.**
  MAPQ is conservative there, so a caller that weighs reads by MAPQ recalls fewer SNPs than with BWA-MEM or minimap2.
- **Many short mapping runs against one large reference.**
  `map` rebuilds the reference index each time, and it uses more memory than minimap2 and much more than BWA-MEM.
- **Any machine that's not Linux on x86-64 (wsl might work but it is untested).**

## Benchmarks and validation

The scripts in `bench/` also need Python 3, hyperfine, Java 17 or newer for Picard, and the libraries the compared programs compile against.

| Platform | Packages |
|---|---|
| Debian, Ubuntu | `sudo apt install python3 hyperfine cmake autoconf automake pkg-config libbz2-dev liblzma-dev libncurses-dev libcurl4-openssl-dev libboost-all-dev libglpk-dev perl openjdk-17-jre-headless samtools bcftools bwa` |
| Fedora, RHEL | `sudo dnf install python3 hyperfine cmake autoconf automake bzip2-devel xz-devel ncurses-devel libcurl-devel boost-devel glpk-devel perl java-21-openjdk-headless samtools bcftools bwa` |
| Arch | `sudo pacman -S --needed python hyperfine cmake bzip2 xz ncurses curl boost glpk perl jre-openjdk-headless` |
| openSUSE | `sudo zypper install python3 hyperfine cmake autoconf automake libbz2-devel xz-devel ncurses-devel libcurl-devel boost-devel glpk-devel perl java-21-openjdk-headless` |
| Nix | `nix-shell -p python3 hyperfine cmake autoconf automake bzip2 xz ncurses curl boost glpk perl jdk21_headless samtools bcftools bwa` |

Before any bench runs, build the compared tools and download the genomes and datasets.

```bash
bench/build_tools.sh                                   # SSHash, SBWT, FMSI, Cuttlefish, Bifrost, KMC, minimap2, BWA and the rest, into temp/
python3 bench/compare_tools_published_sets.py tools    # GGCAT, pbsim2, Winnowmap, mason2 and the other tools of the paper datasets
bench/build_tools.sh scaling                           # hg38 chr21 and chr1 for the scaling bench
bench/build_tools.sh hg38                              # the whole hg38 assembly for the GRCh38 graph benches
python3 bench/compare_tools_published_sets.py list     # each paper dataset, its download size and whether this machine can run it
python3 bench/compare_tools_published_sets.py fetch    # the paper datasets that fit, resumable
```

While `fetch` runs, `python3 bench/compare_tools_published_sets.py status --watch 10` in another terminal shows each dataset's progress, the download rate and the time left.
The E. coli genome some benches use needs no command, because it downloads itself on first use.
With the inputs in place, `run_benchmarks.py` lists and runs the benches.

```bash
python3 bench/run_benchmarks.py list                   # every bench, what it needs and whether it is ready
python3 bench/run_benchmarks.py run --only validations,small_sets,scaling
python3 bench/run_benchmarks.py run --all              # everything, the paper datasets included
```

| Bench | What it does | Time on 16 threads |
|---|---|---|
| `validations` | Equality with the original tools and independent oracles, stress tests, negative controls | 28 min |
| `small_sets` | Head-to-head with SSHash, SBWT, FMSI, Cuttlefish, Bifrost, minimap2 and wfmash on small genomes | 75 min |
| `scaling` | Index, graph and mapping cost from E. coli to human chr1 | 6.5 min to chr21 |
| `published_smoke` | All tools on a slice of each paper dataset, to catch format problems early | 6 min |
| `published` | The datasets and protocols of the compared tools' papers | about a day |
| `dbg`, `wgsim`, `paired`, `mapq`, `real_*`, `micro` | Graph builds on read sets, short-read accuracy, pairing, MAPQ calibration, real reads, hardware counters | varies |

Each bench writes its own section of `bench/results/RESULTS.md` with the date, commit, machine and command line, and a rerun replaces just that section.
`python3 bench/report.py --json` returns the tables as data.
When an input is missing, a bench skips that part or stops and names the command that provides it.
Paper datasets too large for this machine's disk or memory are held back, as `list` shows, and the 1.1 TB `cf2-white-spruce` downloads only with `--only cf2-white-spruce`.
Every script scales its threads and memory cap to the machine, takes `--warmup` and `--runs` for the repetitions, and describes itself with `--help`.
