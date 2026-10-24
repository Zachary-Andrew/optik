"""SAM fields and constants the bench scripts read.

Flag bits and column count from the SAM specification (The SAM/BAM Format Specification Working Group, section 1.4).
HIGH_MAPQ is the confident-placement cut, MAPQ >= 30, that the mapping benches report.
"""

# mandatory columns of an alignment line, QNAME to QUAL
SAM_MIN_COLUMNS = 11

SAM_FLAG_PAIRED = 0x1
SAM_FLAG_PROPER_PAIR = 0x2
SAM_FLAG_UNMAPPED = 0x4
SAM_FLAG_MATE_UNMAPPED = 0x8
SAM_FLAG_REVERSE = 0x10
SAM_FLAG_MATE_REVERSE = 0x20
SAM_FLAG_READ1 = 0x40
SAM_FLAG_READ2 = 0x80
SAM_FLAG_SECONDARY = 0x100
SAM_FLAG_SUPPLEMENTARY = 0x800

HIGH_MAPQ = 30

# optik dbruijn's largest k, since its (k+1)-mer edges must fit 64 bits
DBG_MAX_K = 31

# records that are not a mapped primary alignment
SAM_FLAGS_NOT_PRIMARY = SAM_FLAG_UNMAPPED | SAM_FLAG_SECONDARY | SAM_FLAG_SUPPLEMENTARY

# records that are secondary or supplementary, i.e. not the read's primary line
SAM_FLAGS_SECONDARY_SUPPLEMENTARY = SAM_FLAG_SECONDARY | SAM_FLAG_SUPPLEMENTARY
