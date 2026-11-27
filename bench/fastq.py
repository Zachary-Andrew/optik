"""FASTQ record layout and the two record-level helpers the bench scripts share.

A FASTQ record here is the plain four-line form: @name, sequence, +, quality.
Every read set the benches generate or download uses that form, so counting lines is enough.
"""

FASTQ_LINES_PER_RECORD = 4
FASTQ_SEQ_LINE = 1
FASTQ_QUAL_LINE = 3


def fastq_sequences(lines):
    """The sequence line of every record in `lines`, stripped of its newline."""
    return [line.strip() for line in lines[FASTQ_SEQ_LINE::FASTQ_LINES_PER_RECORD]]


def count_fastq_records(handle):
    """Records left in an open FASTQ text handle, consuming it."""
    return sum(1 for _ in handle) // FASTQ_LINES_PER_RECORD
