"""Describes the record layouts of every format the benches read.

Those are SAM, PAF, GFA, VCF, the truth TSVs and optik's kdbg archive.

Columns are IntEnums, so `fields[SamCol.POS]` keeps working on a plain split line while naming what it reads.
SamRecord goes one step further and wraps a split SAM line with typed accessors.
SAM details follow the SAM/BAM Format Specification Working Group, section 1.4.
"""

import enum

# Mandatory columns of an alignment line, QNAME to QUAL.
SAM_MIN_COLUMNS = 11
# An optional field reads TAG:TYPE:VALUE, so VALUE starts after the 5-character "TG:T:" prefix and ':' follows the 2-character name.
SAM_TAG_PREFIX_LEN = 5
SAM_TAG_NAME_LEN = 2
SAM_MAX_QNAME_LEN = 254

# MAPQ is -10 log10 P(wrong), and SAM_MAX_MAPQ is the largest value optik and minimap2 write.
PHRED_SCALE = 10
SAM_MAX_MAPQ = 60
# The confident-placement cut that the mapping benches report.
HIGH_MAPQ = 30

# optik dbruijn's largest k, since its (k+1)-mer edges must fit 64 bits.
DBG_MAX_K = 31


class SamFlag(enum.IntFlag):
    """FLAG bits of a SAM record."""

    PAIRED = 0x1
    PROPER_PAIR = 0x2
    UNMAPPED = 0x4
    MATE_UNMAPPED = 0x8
    REVERSE = 0x10
    MATE_REVERSE = 0x20
    READ1 = 0x40
    READ2 = 0x80
    SECONDARY = 0x100
    SUPPLEMENTARY = 0x800
    # Records that are not a mapped primary alignment.
    NOT_PRIMARY = UNMAPPED | SECONDARY | SUPPLEMENTARY
    # Secondary or supplementary, meaning some other line is the read's primary.
    SECONDARY_SUPPLEMENTARY = SECONDARY | SUPPLEMENTARY


class SamCol(enum.IntEnum):
    """0-based columns of a SAM alignment line, with optional tags starting at TAGS."""

    QNAME = 0
    FLAG = 1
    RNAME = 2
    POS = 3
    MAPQ = 4
    CIGAR = 5
    RNEXT = 6
    PNEXT = 7
    TLEN = 8
    SEQ = 9
    QUAL = 10
    TAGS = SAM_MIN_COLUMNS


class SamRecord(tuple):
    """A SAM alignment line split on tabs.

    It is still a tuple, so existing `fields[SamCol.X]` code works unchanged, but the accessors say what they read.
    parse() splits exactly like line.split("\\t"), so a trailing newline stays on the last field.
    """

    __slots__ = ()

    @classmethod
    def parse(cls, line):
        return cls(line.split("\t"))

    qname = property(lambda self: self[SamCol.QNAME])
    rname = property(lambda self: self[SamCol.RNAME])
    cigar = property(lambda self: self[SamCol.CIGAR])
    rnext = property(lambda self: self[SamCol.RNEXT])
    seq = property(lambda self: self[SamCol.SEQ])
    qual = property(lambda self: self[SamCol.QUAL])
    flag = property(lambda self: int(self[SamCol.FLAG]))
    pos = property(lambda self: int(self[SamCol.POS]))
    mapq = property(lambda self: int(self[SamCol.MAPQ]))
    pnext = property(lambda self: int(self[SamCol.PNEXT]))
    tlen = property(lambda self: int(self[SamCol.TLEN]))
    tags = property(lambda self: self[SamCol.TAGS:])

    @property
    def is_unmapped(self):
        return bool(self.flag & SamFlag.UNMAPPED)

    @property
    def is_reverse(self):
        return bool(self.flag & SamFlag.REVERSE)

    @property
    def is_primary(self):
        """Mapped, and neither secondary nor supplementary."""
        return not self.flag & SamFlag.NOT_PRIMARY


class PafCol(enum.IntEnum):
    """0-based columns of a PAF line (minimap2 man page)."""

    QNAME = 0
    QLEN = 1
    QSTART = 2
    QEND = 3
    STRAND = 4
    TNAME = 5
    TLEN = 6
    TSTART = 7
    TEND = 8
    MATCHES = 9
    BLOCK_LEN = 10
    MAPQ = 11
    TAGS = 12


class GfaSegmentCol(enum.IntEnum):
    """0-based columns of a GFA S line, after the record type."""

    NAME = 1
    SEQ = 2
    TAGS = 3


class GfaLinkCol(enum.IntEnum):
    """0-based columns of a GFA L line, after the record type."""

    FROM = 1
    FROM_ORIENT = 2
    TO = 3
    TO_ORIENT = 4
    OVERLAP = 5


class VcfCol(enum.IntEnum):
    """0-based columns of a VCF data line."""

    CHROM = 0
    POS = 1
    REF = 3
    ALT = 4
    QUAL = 5
    INFO = 7


class TruthCol(enum.IntEnum):
    """Columns of gen_reads.py's ground-truth TSV: read name, reference, 0-based position, strand, read length."""

    NAME = 0
    REF = 1
    POS = 2
    STRAND = 3
    LEN = 4


# optik dbruijn -s node table (OKKDBG), little-endian:
# magic | u8 k | u64 node count n | n x u64 sorted canonical k-mers | n x u32 abundances | n x u8 edge masks (in << KDBG_IN_MASK_SHIFT | out)
KDBG_MAGIC = b"OKKDBG"
KDBG_K_OFFSET = len(KDBG_MAGIC)
KDBG_COUNT_OFFSET = KDBG_K_OFFSET + 1
KDBG_HEADER_BYTES = KDBG_COUNT_OFFSET + 8
KDBG_KEY_BYTES = 8
KDBG_ABUNDANCE_BYTES = 4
KDBG_MASK_BYTES = 1
KDBG_NODE_BYTES = KDBG_KEY_BYTES + KDBG_ABUNDANCE_BYTES + KDBG_MASK_BYTES
KDBG_IN_MASK_SHIFT = 4
