"""Dataset catalogue of compare_tools_published_sets.py.

Each dataset the published-sets bench can run is registered with ds(), which records what it is, the files to fetch, how its paper used it, the published counts and memory figures that serve as gates, and the paper it comes from.
The module also holds the paths of the bench's working directories and the pinned helper-tool versions.
"""

import os

from external_oracles import BWA_BIN

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


T = os.path.join(ROOT, "temp")

TOOLS = os.path.join(T, "published_sets_tools")

# Helper tools are taken at their latest stable version, meaning default-branch HEAD or the latest release where the tool ships that way, as k8's prebuilt binary and SeqAn 2's mason2 do.
# The papers fix what is measured and with which parameters, not these versions.
K8_URL = "https://github.com/attractivechaos/k8/releases/download/v1.2/k8-1.2.tar.bz2"

BIWFA_REV = "6289fd23f2a908097b4532ec67fad2f1188c7a3e"

SEQAN_TAG = "seqan-v2.5.3"  # latest SeqAn 2 release, which provides mason2

HG38 = os.path.join(
    T, "scaling/hg38_primary.fa"
)  # UCSC's 25 primary chromosomes of GRCh38

BWA_IDX = os.path.join(T, "real/bwaidx/hg38.fa")  # bwa index of the same file (symlink)

BWA = BWA_BIN or os.path.join(T, "bwa/bwa")  # from PATH, else bench/build_tools.sh bwa

Z_SSHASH = "17582116"  # "Eulertigs for benchmarking kmer dictionaries" (SSHash 2026 data availability)

Z_FMSI = "14722244"  # "Genomic datasets used for evaluation of k-mer representations and indexes"

Z_SBWT = "6577997"  # "3682 E. coli assemblies from NCBI" (SBWT repository, ecoli_dataset/README.md)

SSHASH26 = "Pibiri & Patro 2026, Bioinformatics btag264 (SSHash v6)"

FMSI25 = "Sladky, Vesely & Brinda 2025, Bioinformatics Advances vbaf290 (FMSI)"

SBWT23 = (
    "Alanko, Puglisi & Vuohtoniemi 2023, ACDA (SBWT); bioRxiv 10.1101/2022.05.19.492613"
)

CF2 = "Khan, Kokot, Deorowicz & Patro 2022, Genome Biology 23:190 (Cuttlefish 2)"

CF3 = "Khan et al. 2026, bioRxiv 10.1101/2025.02.02.636161 (Cuttlefish 3)"

BIFROST = "Holley & Melsted 2020, Genome Biology 21:249 (Bifrost)"

BWAMEM = "Li 2013, arXiv:1303.3997 (BWA-MEM)"

MM2_18 = "Li 2018, Bioinformatics 34:3094 (minimap2)"

MM2_21 = "Li 2021, Bioinformatics 37:4572, btab705 (minimap2 2.22 strategies)"

BIWFA = "Marco-Sola et al. 2023, Bioinformatics 39:btad074 (BiWFA)"

KMC3 = "Kokot, Dlugosz & Deorowicz 2017, Bioinformatics 33:2759 (KMC 3)"


def zen(rec, key):
    """File reference to `key` in Zenodo record `rec`, as used in a dataset's `files`."""
    return {"zenodo": rec, "key": key}


# Dataset manifest
# Maps an id to its dataset.
# `title` says what the dataset is and `files` lists what is fetched.
# `use` describes how the paper used it and how the bench uses it.
# `kmers31` is the published distinct 31-mer count, and `mem` maps a tool to the (GB, source) memory figure a paper reports for this dataset.
D = {}


def ds(id_, **kw):
    """Register dataset id_ in D, defaulting to no files and tier M.

    `tscale` is the first thread count of a thread-scaling sweep.
    The bench runs hw.thread_points() from there up to its --threads, by default every usable CPU, so the sweep grows with the machine."""
    kw.setdefault("files", [])
    kw.setdefault("tier", "M")
    D[id_] = kw


# Query datasets, whole genomes and pangenomes, from the SSHash 2026 collections
_SS_Q = {  # collection -> (title, high-hit query file, SSHash -m at k=31, distinct 31-mers, tier)
    "cod": ("Cod (Gadus morhua) genome", "SRR12858649.fastq.gz", 20, 502_465_200, "M"),
    "kestrel": (
        "Kestrel (Falco tinnunculus) genome",
        "SRR11449743_1.fastq.gz",
        20,
        1_150_399_205,
        "M",
    ),
    "human": (
        "Human genome (GRCh38, Ensembl toplevel)",
        "SRR5833294.fastq.gz",
        21,
        2_505_678_680,
        "L",
    ),
    "ncbi-virus": (
        "18,836 RefSeq virus genomes (pangenome)",
        "ncbi-queries.fastq.gz",
        19,
        376_205_185,
        "M",
    ),
    "se": (
        "534,751 Salmonella enterica genomes (pangenome)",
        "SRR27871075_1.fastq.gz",
        21,
        894_310_084,
        "L",
    ),
}

for c, (title, q, m, n, tier) in _SS_Q.items():
    ds(
        f"sshash26-{c}",
        section="query",
        paper=SSHASH26,
        tier=tier,
        title=title + " eulertigs, k = 31",
        source=f"zenodo.org/records/{Z_SSHASH} (paper: 'Data availability'); query pairing from "
        "kmer_sets_benchmark/script/query_sshash_indexes.sh",
        files=[zen(Z_SSHASH, f"{c}.k31.eulertigs.fa.gz"), zen(Z_SSHASH, q)],
        input=f"{c}.k31.eulertigs.fa.gz",
        stream=q,
        k=[31],
        sshash_m=m,
        kmers31=n,
        use="index the eulertigs; 1e6 positive (half RC) and 1e6 random k-mers; stream the high-hit FASTQ",
    )

ds(
    "sshash26-hprc",
    section="query",
    paper=SSHASH26,
    tier="L",
    title="Human pangenome (HPRC, 472 haplotypes) eulertigs, k = 31",
    source=f"zenodo.org/records/{Z_SSHASH}; HPRC human472.agc from zenodo.org/records/14854401 (paper)",
    files=[
        zen(Z_SSHASH, "hprc.k31.eulertigs.fa.gz"),
        zen(Z_SSHASH, "SRR5833294.fastq.gz"),
    ],
    input="hprc.k31.eulertigs.fa.gz",
    stream="SRR5833294.fastq.gz",
    k=[31],
    sshash_m=21,
    kmers31=3_718_120_949,
    mem={
        "fmsi": (
            75.2,
            (
                "SSHash 2026 Results: FMSI 'required 70 and 125 GiB of RAM on HPRC for k = 31 and "
                "k = 63' (70 GiB = 75.2 GB)"
            ),
        )
    },
    use="as the other collections (3.7 B distinct 31-mers)",
)

ds(
    "sweep-chr1",
    section="query",
    paper=SSHASH26,
    tier="M",
    title="Human chromosome 1 (GRCh38), the matched k x threads sweep",
    source="chromosome 1 of the GRCh38 genome SSHash 2026 indexes as its 'human' collection (UCSC copy, "
    "bench/build_tools.sh scaling); queried with the same deposit's high-hit FASTQ SRR5833294",
    local=os.path.join(T, "scaling/chr1.fa"),
    stream="SRR5833294.fastq.gz",
    stream_from="sshash26-human",
    k=[15, 19, 21, 23, 27, 31],
    spss=True,
    kscale=[15, 23, 31],
    tscale=1,
    use="SPSS by GGCAT at every k, then every index at the same k and thread count; the comparison set "
    "of the sweep (the deposit's eulertigs exist only at k = 31)",
)

ds(
    "sshash26-k63",
    section="query",
    paper=SSHASH26,
    tier="X",
    title="The SSHash 2026 sets at k = 63",
    excluded="the paper's k = 63 sets: optik supports k <= 32, so there is nothing to compare",
    source=f"zenodo.org/records/{Z_SSHASH}",
)

# Query datasets from the FMSI deposit (unitigs at k = 128, recomputed as an SPSS at each k)
_FM = {  # id -> (zenodo key, title, paper's name, tier, 10%-subsampled k = 23 file)
    "spneumo": (
        "spneumo_pangenome_RASE_db.unitigs_k128.fa.xz",
        "616 S. pneumoniae genomes (pangenome)",
        "sp-pg (Croucher et al. 2015)",
        "S",
        "spneumo_pangenome_RASE_db.unitigs_k128_subsampled_k23_r0.1.fa.xz",
    ),
    "sars-cov-2": (
        "sars-cov-2_pangenome_gisaid.unitigs_k128.fa.xz",
        "SARS-CoV-2 GISAID genomes (pangenome)",
        "sc2-pg (GISAID; the derived unitigs are deposited)",
        "S",
        "sars-cov-2_pangenome_gisaid.unitigs_k128_subsampled_k23_r0.1.fa.xz",
    ),
    "ecoli-hq": (
        "ecoli_pangenome_661k_HQ.unitigs_k128.fa.xz",
        "E. coli high-quality genomes of the 661k collection",
        "ec-pg-hq (Blackwell 2021)",
        "M",
        "ecoli_pangenome_661k_HQ.unitigs_k128_subsampled_k23_r0.1.fa.xz",
    ),
    "microbiome": (
        "human_microbiome_illumina_srs063932.unitigs_minfreq2_k32.fa.xz",
        "Human microbiome Illumina reads (k-mers seen twice)",
        "mtg-ilm (SRX023459)",
        "M",
        "human_microbiome_illumina_srs063932.unitigs_minfreq2_k32_subsampled_k23_r0.1.fa.xz",
    ),
    "ecoli-all": (
        "ecoli_pangenome_661k_all.unitigs_k128.fa.xz",
        "All E. coli genomes of the 661k collection",
        "ec-pg-all",
        "L",
        "ecoli_pangenome_661k_all.unitigs_k128_subsampled_k23_r0.1.fa.xz",
    ),
    "rnaseq": (
        "human_rnaseq_srx348811.unitigs_minfreq2_k32.fa.xz",
        "Human RNA-seq reads (k-mers seen twice)",
        "rna-ilm (SRX348811)",
        "L",
        "human_rnaseq_srx348811.unitigs_minfreq2_k32_subsampled_k23_r0.1.fa.xz",
    ),
    "chm13": (
        "human_genome_assembly_chm13.v2.0.fa.xz",
        "Human T2T-CHM13 v2.0 assembly",
        "hg-t2t (Nurk 2022)",
        "L",
        "human_genome_assembly_chm13.v2.0_subsampled_k23_r0.1.fa.xz",
    ),
    "hg-ilm": (
        "human_genome_illumina.unitigs_minfreq2_k32.fa.xz",
        "Human genome Illumina reads (k-mers seen twice)",
        "hg-ilm (SRX016231)",
        "L",
        None,
    ),
}

_FM_UNITIGS = {  # id -> why GGCAT's unitigs stand in for its eulertigs as the SPSS
    "ecoli-all": "eulertigs: GGCAT issue #65, open; its eulertig joiner aborts on this collection "
    "(a wrapped length in sequences_joiner.rs), 7 of 8 runs at k = 23 and 31, 1 or 8 threads",
}

_FM_EXCLUDED = {  # id -> why the deposit cannot be used
    "ecoli-hq": "the deposited ecoli_pangenome_661k_HQ.unitigs_k128.fa.xz is a truncated xz stream: its size "
    "(528,744,448 B) and MD5 (869f308dc8b07085837f86ccfd8f2369) are those Zenodo lists, but it "
    "ends inside unitig 30,698,059 (180 of LN:i:242 bases) with no end-of-stream marker, so the "
    "full set cannot be recovered (checked 2026-10-17)",
}

for c, (key, title, name, tier, sub) in _FM.items():
    ds(
        f"fmsi25-{c}",
        section="query",
        paper=FMSI25,
        tier=tier,
        title=title + ", k = 23 and 31",
        **({"excluded": _FM_EXCLUDED[c]} if c in _FM_EXCLUDED else {}),
        **({"spss_unitigs": _FM_UNITIGS[c]} if c in _FM_UNITIGS else {}),
        source=f"doi.org/10.5281/zenodo.{Z_FMSI} (paper: 'The datasets used in the experiments are deposited "
        "at Zenodo'); paper's name: " + name,
        files=[zen(Z_FMSI, key)] + ([zen(Z_FMSI, sub)] if sub else []),
        input=key,
        subsampled=sub,
        k=[23, 31],
        spss=True,
        sshash_m_rule="fmsi",
        # FMSI's memory k-sweep on its large sets uses k = 15, 23 and 31 (Table 4, Fig. 2), which become the plan's k points.
        **({"kscale": [15, 23, 31]} if c == "ecoli-all" else {}),
        use="SPSS at each k (the paper's GGCAT step), then every index; the 10%-subsampled k = 23 set "
        "indexed as deposited; FMSI's query sets are in its supplement (not public text), so the "
        "SSHash 2026 positive/negative protocol is used",
    )

ds(
    "sbwt23-ecoli3682",
    section="query",
    paper=SBWT23,
    tier="M",
    title="3,682 E. coli genomes (pangenome), k = 31",
    source=f"zenodo.org/record/{Z_SBWT} (SBWT repo ecoli_dataset/README.md; paper Sec. 6.1: 3682 E. coli "
    "genomes, 257,785,486 distinct 31-mers)",
    files=[zen(Z_SBWT, "coli3682_dataset.tar.gz")],
    k=[31],
    spss=True,
    kmers31=257_785_486,
    # The paper counts the genomes with each non-ACGT base deleted inside its record (pieces joined, not split) and k-mers as written (forward).
    # KMC on 2026-10-18 gave exactly 257,785,486 forward k-mers this way.
    # Splitting at non-ACGT gives 253,718,096 forward and 170,648,610 canonical.
    delete_nonacgt=True,
    kmers31_forward=True,
    use="non-ACGT deleted within each record (paper), SPSS at k = 31, then every index; 1e6 positive "
    "and random k-mers",
)

ds(
    "sbwt23-gut",
    section="query",
    paper=SBWT23,
    tier="L",
    title="Human gut metagenome reads ERR5035349, k = 31",
    source="SBWT paper Sec. 6.1: '17 million Illumina HiSeq 2500 sequence reads of length 502 sampled from the "
    "human gut (SRA identifier ERR5035349)', 2,771,108,020 distinct 31-mers",
    files=[{"ena": "ERR5035349"}],
    reads=True,
    k=[31],
    spss=True,
    kmers31=2_771_108_020,
    # The paper's "reads of length 502" are the run's 17.3 M spots with mates joined (2 x 251, SRA's spot), and its count is of k-mers as written (forward), not canonical.
    # KMC on ENA's FASTQ (2026-10-16) found 2,761,523,935 forward k-mers in the joined spots (99.65 %).
    # With mates apart it found 2,258,199,229 forward and 1,873,574,701 canonical.
    # Non-ACGT handling does not explain the remaining 0.35 % (only 3.98 M windows touch such a base), so the gate allows 0.5 %.
    join_mates=True,
    kmers31_forward=True,
    kmers31_tol=0.005,
    use="every k-mer of the reads (no abundance cutoff, as the paper), SPSS, every index",
)

ds(
    "sbwt23-sarscov2",
    section="query",
    paper=SBWT23,
    tier="X",
    title="14.6 M SARS-CoV-2 genomes from NCBI",
    excluded="14,641,164 SARS-CoV-2 genomes 'downloaded from NCBI datasets' with no accession list or date",
    source="paper Sec. 6.1",
)

# De Bruijn graph datasets of Cuttlefish 2
ds(
    "cf2-celegans-ref",
    section="dbg",
    paper=CF2,
    tier="S",
    title="C. elegans VC2010 genome, k = 27",
    source="CF2 'roundworm reference (Caenorhabditis elegans, VC2010) [68]'; [68] Yoshimura et al. 2019 "
    "Genome Res: assembly deposited in ENA PRJEB28388 = WGS set UNSB01",
    files=[
        {
            "url": "https://ftp.ebi.ac.uk/pub/databases/ena/wgs/public/uns/UNSB01.fasta.gz"
        }
    ],
    input="UNSB01.fasta.gz",
    k=[27],
    cutoff=1,
    kind="ref",
    tscale=1,
    use="compacted graph of the reference, f0 = 1 (CF2 Table 3); thread scaling from 1 thread (a small "
    "genome: the single-threaded point is affordable here, not on HG004)",
)

ds(
    "cf2-celegans-reads",
    section="dbg",
    paper=CF2,
    tier="M",
    title="C. elegans Illumina reads DRR008444, k = 27",
    source="CF2: 'a roundworm read set (ENA DRR008444) ... 5.6 GB compressed FASTQ'; f0 = 12 at k = 27",
    files=[{"ena": "DRR008444"}],
    k=[27],
    cutoff=12,
    kind="reads",
    tscale=1,
    use="compacted graph of both mates, abundance >= 12 (CF2 Table 3); thread scaling from 1 thread",
)

ds(
    "cf2-grch38",
    section="dbg",
    paper=CF2,
    tier="L",
    title="Human genome GRCh38, k = 27",
    kscale=[15, 21, 27, 31],
    source="CF2: 'a human reference (Homo sapiens, GRCh38)'; local UCSC GRCh38 primary chromosomes",
    local=HG38,
    k=[27],
    cutoff=1,
    kind="ref",
    use="compacted graph of GRCh38, f0 = 1",
)

# Read sets of Cuttlefish 2 Table 1, with f0 at k = 27 taken from its footnote ("human: 14 (k = 27) ..., human RNA-seq, gut microbiome and soil: 2, white spruce: 11 (k = 27)").
# The paper's sizes are in GiB, and ENA's byte counts agree with them, 158 GB = 147 GiB for HG004, except for the soil set, whose note explains why.
# Table 1's peak memory figures in GB, at k = 27 with 8 threads for human and 16 otherwise, are for Bifrost and Cuttlefish 2's default.
_CF2_READS = {
    "hg004": (
        "Human HG004 Illumina reads (70x)",
        [f"SRR{n}" for n in range(3440461, 3440496)],
        14,
        (
            "'a human read set (NIST HG004) ... SRA3440461-95, 148 GB compressed FASTQ' = SRR3440461-SRR3440495 "
            "(ENA sample 'NIST HG004 NA24143')"
        ),
        None,
        {"bifrost": 48.5, "cuttlefish": 3.2},
    ),
    "rnaseq-1kg": (
        "1000 Genomes small-RNA-seq reads (465 samples)",
        ["ERP001941"],
        2,
        (
            "'RNA sequencing dataset (ENA PRJEB3365) of 465 human lymphoblastoid cell line samples ... "
            "ERP001941, 140 GB compressed FASTQ' (ENA: 492 runs)"
        ),
        None,
        {"bifrost": 7.2, "cuttlefish": 2.9},
    ),
    "gut": (
        "Human gut microbiome reads (9 individuals)",
        ["ERP115863"],
        2,
        (
            "'a gut microbiome read set (ENA PRJEB33098) from nine individuals ... ERP115863, 45 GB compressed "
            "FASTQ'"
        ),
        None,
        {"bifrost": 38.1, "cuttlefish": 3.5},
    ),
    "soil": (
        "Iowa corn soil metagenome reads",
        ["SRX100357", "SRX099904", "SRX099905", "SRX099906"],
        2,
        "'soil metagenome read set (Iowa Corn) ... SRX100357 and SRX099904-06, 152 GB compressed FASTQ'",
        (
            "ENA serves these four experiments as 87.7 GB of FASTQ.gz, not the paper's 152 GB; the "
            "accessions are the paper's, the difference is recorded, not reconciled"
        ),
        {"bifrost": 274.1, "cuttlefish": 19.2},
    ),
    "white-spruce": (
        "White spruce (Picea glauca) reads",
        ["PRJNA83435"],
        11,
        (
            "'a white spruce read set (NCBI PRJNA83435) ... SRA056234, 1.14 TB compressed FASTQ' "
            "(ENA: 61 runs, 1,147 GB)"
        ),
        None,
        {"cuttlefish": 14.0},
    ),
}

# Why white spruce runs only when named, shown by `list` and `status`
OPT_IN_WHITE_SPRUCE = (
    "opt-in: 1,147 GB of reads, fetched and run only when named, as in "
    "`fetch --only cf2-white-spruce`"
)

for c, (title, accs, f0, src, note, mem) in _CF2_READS.items():
    ds(
        f"cf2-{c}",
        section="dbg",
        paper=CF2,
        tier="L",
        title=title + ", k = 27",
        source="CF2 Table 1: " + src,
        files=[{"ena": a} for a in accs],
        k=[27],
        cutoff=f0,
        kind="reads",
        note=note,
        mem={t: (gb, "CF2 Table 1, peak at k = 27") for t, gb in mem.items()},
        skip_tools={
            "bifrost": "CF2 did not run Bifrost here: it cannot take f0 = 11 (fixed f0 = 2)"
        }
        if c == "white-spruce"
        else {},
        use=f"compacted graph of all runs, f0 = {f0} at k = 27 (CF2 Table 1)",
        # 1.1 TB of reads is never fetched or run by default, whatever the free disk.
        **({"opt_in": OPT_IN_WHITE_SPRUCE} if c == "white-spruce" else {}),
    )

# Cuttlefish 2's scaling experiments use HG004 "downsampled ... from 70x to 20x coverage".
# Thread scaling runs at k = 27, f0 = 4 over 1-32 threads (Fig. 2), and k scaling at f0 = 5 with 8 threads over k = 27..125 (Table 4).
# The paper does not identify its subsample, so this bench takes the prefix of whole runs of SRR3440461-95 that comes nearest 20/70 of the set's compressed bytes.
# Run sizes vary.
# On 2026-10-20 that was the first 28 runs, 32.6% of the bytes and about 23x, and the bench prints the share it uses.
# The thread points run up to the machine's usable CPUs, and the k points are limited to optik's k <= 31, with odd k for Cuttlefish.
ds(
    "scale-dbg-hg004-20x",
    section="dbg",
    paper=CF2 + "; " + CF3,
    tier="L",
    as_="cf2-hg004",
    files=[],
    title="HG004 reads at 20x, thread scaling (k = 27, f0 = 4)",
    source="CF2 Fig. 2 protocol (HG004 downsampled 70x -> 20x, k = 27, f0 = 4); runs: the prefix of whole "
    "runs of SRR3440461-95 nearest 20/70 of the compressed bytes (28 runs, 32.6%, ~23x)",
    k=[27],
    cutoff=4,
    kind="reads",
    runs_frac=20 / 70,
    no_c2_round=True,
    tscale=4,
    tools=["optik", "cuttlefish"],
    use="optik and Cuttlefish from 4 threads up to every usable CPU, same k, cutoff and memory budget (CF2 Fig. 2 starts "
    "at 1 thread; on a human read set that point takes many hours per tool, so the single-threaded "
    "scaling is measured on the C. elegans sets)",
)

ds(
    "kscale-dbg-hg004-20x",
    section="dbg",
    paper=CF2 + "; " + CF3,
    tier="L",
    as_="cf2-hg004",
    files=[],
    title="HG004 reads at 20x, k scaling (f0 = 5, 8 threads)",
    source="CF2 Table 4 protocol (HG004 20x, f0 = 5, 8 threads), k limited to optik's k <= 31",
    k=[27],
    cutoff=5,
    kind="reads",
    runs_frac=20 / 70,
    no_c2_round=True,
    kscale=[15, 21, 27, 31],
    tools=["optik", "cuttlefish"],
    use="optik and Cuttlefish at k = 15, 21, 27, 31 at 8 threads",
)

ds(
    "cf2-humgut",
    section="dbg",
    paper=CF2 + "; " + CF3,
    tier="X",
    title="HumGut 30,691 gut genomes",
    excluded="HumGut (30,691 genomes, Hiseni et al. 2021): its URL (arken.nmbu.no/~larssn/humgut/) now "
    "serves only HumGut2, a different collection",
    source="CF2 Table 2; CF3 Sec. 5.1",
)

ds(
    "cf3-bacteria",
    section="dbg",
    paper=CF3,
    tier="X",
    title="309K Salmonella / 661K bacteria (colored)",
    excluded="colored graphs (309K Salmonella, 661K bacteria): optik builds uncolored graphs; the 10,000-"
    "bacteria scaling subset is not identified in the paper",
    source="CF3 Sec. 5.1-5.3",
)

ds(
    "bifrost-na12878",
    section="dbg",
    paper=BIFROST,
    tier="X",
    title="NA12878 reads subsampled to 30x",
    excluded="NA12878 GIAB reads 'downsampled from 300-fold to 30-fold': the paper names neither the "
    "runs nor the subsampling, and its Salmonella set is a colored graph",
    source="Bifrost Results",
)

ds(
    "kmc3-sets",
    section="count",
    paper=KMC3,
    tier="X",
    title="KMC 3 counting sets",
    excluded="KMC 3's datasets are described only in its supplementary material (no accessions in the "
    "text); KMC is an oracle here, not a compared tool",
    source="KMC 3 Results",
)

# Mapping datasets
ds(
    "bwamem-wgsim-hg38",
    section="map",
    paper=BWAMEM,
    tier="L",
    title="1 M simulated 101 bp read pairs from the human genome (wgsim)",
    source="BWA-MEM Fig. 1: '10^6 pairs of 101bp reads are simulated from the human reference genome "
    "using wgsim with 1.5% substitution errors and 0.2% indel variants. The insert size follows "
    "N(500, 50^2)'; correct = within 20 bp; MAPQ <= 3 excluded",
    local=HG38,
    sim="wgsim",
    pairs=1_000_000,
    length=101,
    use="wgsim -N 1e6 -1 101 -2 101 -d 500 -s 50 -e 0.015 -r 0.002 -R 1.0 (all variants indels); "
    "mapped paired and single-end; wgsim_eval.pl at 20 bp",
)

ds(
    "mm2-pbsim2-chm13",
    section="map",
    paper=MM2_21,
    tier="L",
    title="152,713 simulated Nanopore reads from T2T-CHM13 v1.1 (pbsim2)",
    source="btab705: '152 713 reads were simulated from the CHM13 telomere-to-telomere assembly v1.1 "
    "(AC: GCA_009914755.3) with pbsim2 --hmm_model R94.model --length-min 5000 --length-mean 20000 "
    "--accuracy-mean 0.95'; 'Alignments of mapping quality 10 or higher were evaluated by "
    "paftools.js mapeval'",
    files=[{"ncbi": "GCA_009914755.3"}],
    sim="pbsim2",
    needs=["pbsim2/src/pbsim", "k8"],
    use="pbsim2 at depth 1 (152,713 x 20 kb ~ 1x of 3.05 Gbp; the paper gives the count, not the depth), "
    "names rewritten by paftools.js pbsim2fq, mapped by optik and minimap2 -x map-ont, mapeval",
)

ds(
    "bwamem-ecoli-asm",
    section="map",
    paper=BWAMEM,
    tier="S",
    title="E. coli MG1655 genome aligned to strain 536",
    source="BWA-MEM Results: 'we aligned the E. coli K-12 MG1655 substrain (AC:NC_000913) against the 536 "
    "strain (AC:NC_008253)'",
    files=[{"efetch": "NC_000913.3"}, {"efetch": "NC_008253.1"}],
    use="MG1655 as the query, 536 as the reference; substitutions counted by replaying each tool's CIGAR "
    "(the paper: nucmer 105,505, BWA-MEM 104,321)",
)

ds(
    "mm2-mason2",
    section="map",
    paper=MM2_18,
    tier="L",
    title="10 M simulated 150 bp read pairs from GRCh38 (mason2)",
    source="minimap2 2018 Fig. 1b: '10 million pairs of 150 bp reads were simulated using mason2 with option "
    "--illumina-prob-mismatch-scale 2.5' from the GRCh38 primary assembly; 'Short-read aligners were "
    "run under the default setting except for changing the maximum fragment length to 800 bp'",
    local=HG38,
    needs=["mason_simulator", "k8"],
    pairs=10_000_000,
    use="mason_simulator -n 1e7 --illumina-read-length 150 --illumina-prob-mismatch-scale 2.5 with its SAM "
    "truth, paftools.js mason2fq, optik / minimap2 -ax sr / bwa mem paired, paftools.js mapeval",
    note="the paper's 800 bp maximum fragment length is minimap2's default (-F 800); bwa mem can set it only "
    "through -I, which also fixes the insert mean and sd (the simulator's truth), so bwa and optik "
    "(no such option) infer the insert size from the reads, their defaults",
)

ds(
    "sweep-pbsim2-chm13",
    section="map",
    paper=MM2_21,
    tier="M",
    as_="mm2-pbsim2-chm13",
    subset=20_000,
    title="20,000 of the pbsim2 CHM13 v1.1 reads, the matched k x threads sweep",
    source="every n-th read of mm2-pbsim2-chm13 (btab705's simulation, seed 7), 20,000 over the whole genome",
    files=[],
    needs=["pbsim2/src/pbsim", "k8"],
    k=[15, 19, 21, 23, 27],
    kscale=[15, 21, 27],
    tscale=1,
    use="optik and minimap2 -ax map-ont at the same k (w = 10) and thread count; paftools.js mapeval",
)

ds(
    "sweep-mason2-hg38",
    section="map",
    paper=MM2_18,
    tier="M",
    as_="mm2-mason2",
    pairs=200_000,
    title="200,000 mason2 read pairs from GRCh38, the matched k x threads sweep",
    source="minimap2 2018's mason2 simulation (--illumina-prob-mismatch-scale 2.5, 150 bp pairs) at 200,000 pairs",
    local=HG38,
    needs=["mason_simulator", "k8"],
    k=[15, 19, 21, 23, 27],
    kscale=[15, 21, 27],
    tscale=1,
    use="optik and minimap2 -ax sr at the same k (w = 10) and thread count, bwa mem (no comparable k) at the "
    "same thread count; paftools.js mapeval",
)

ds(
    "mm2-winnowmap-truth",
    section="map",
    paper=MM2_21,
    tier="L",
    title="1.39 M CHM13 HiFi reads judged against Winnowmap2",
    source="btab705 [winno-cmp]: '1.39 million CHM13 HiFi reads from SRR11292121 were mapped against the same "
    "CHM13 assembly. 99.3% of them were mapped by Winnowmap2 at mapping quality 10 or higher and were "
    "taken as ground truth to evaluate minimap2 and lra with paftools.js pafcmp'",
    files=[{"ena": "SRR11292121"}, {"ncbi": "GCA_009914755.3"}],
    needs=["winnowmap", "meryl", "k8"],
    use="Winnowmap2 (its README's repetitive-k-mer recipe, k = 15, -x map-pb) as the truth at MAPQ >= 10; "
    "optik and minimap2 -x map-hifi compared by paftools.js pafcmp",
)

ds(
    "mm2-real-short",
    section="map",
    paper=MM2_18,
    tier="X",
    title="CHM1/CHM13 Illumina reads ERR1341796",
    excluded="ERR1341796 (95 GB) is judged by 'an independent truth variant dataset' built from the CHM1 "
    "and CHM13 assemblies (paper Sec. 3.3), which the paper's text does not name or deposit",
    source="minimap2 2018 Sec. 3.3",
)

ds(
    "mm2-pbsim-chm1",
    section="map",
    paper=MM2_18,
    tier="X",
    title="pbsim SMRT reads with a CHM1 error profile",
    excluded="the pbsim error profile (http://bit.ly/chm1p5c3) redirects to datasets.pacb.com, whose host "
    "name no longer resolves (checked 2026-10-12)",
    source="minimap2 2018 Fig. 1a",
)

# Pairwise alignment datasets
ds(
    "biwfa-ont-ul48",
    section="align",
    paper=BIWFA,
    tier="S",
    title="48 Nanopore read pairs of 500 kbp to 1 Mbp",
    source="BiWFA Sec. 3.1: 'ONT MinION reads from Bowden et al. (2019) ... at least 500 kbp long, for a "
    f"total number of 48 sequence pairs'; file from github.com/smarco/BiWFA-paper @ {BIWFA_REV[:7]}",
    files=[
        {
            "url": f"https://github.com/smarco/BiWFA-paper/raw/{BIWFA_REV}/evaluation/data/"
            "ONT_MinION_UL.500kbps.48.tar.gz"
        }
    ],
    use="all 48 pairs, optik's BiWFA vs WFA2-lib's (MemoryUltralow), penalties 4,6,2",
)

ds(
    "biwfa-sim",
    section="align",
    paper=BIWFA,
    tier="S",
    title="Simulated pairs of 100 kbp to 2 Mbp, 10% and 20% error",
    source="BiWFA Table 1 (100 kbp, 1 Mbp, 2 Mbp at 10% and 20%); generator = WFA2-lib generate_dataset "
    "(BiWFA-paper evaluation/2_Simulations.md: generate_dataset -l LEN -e ERR -n 1)",
    sim="generate_dataset",
    grid=[
        (100000, 0.10),
        (100000, 0.20),
        (1000000, 0.10),
        (1000000, 0.20),
        (2000000, 0.10),
        (2000000, 0.20),
    ],
    use="one pair per cell as in the paper's script; optik's BiWFA vs WFA2-lib's (MemoryUltralow)",
)

ds(
    "biwfa-hprc-chr12",
    section="align",
    paper=BIWFA,
    tier="X",
    title="1,312 HG002 chr12 Nanopore read pairs",
    excluded="the 1312 HG002 chr12 PromethION pairs were prepared on the authors' cluster and are not "
    "deposited",
    source="BiWFA Sec. 3.1; BiWFA-paper evaluation/1_RealData.md",
)

TIERS = {"S": 0, "M": 1, "L": 2, "X": 9}

QUERY_TOOLS = ("optik-spss", "optik-fm", "sshash", "sbwt", "fmsi")

DBG_TOOLS = ("optik", "cuttlefish", "bifrost")
