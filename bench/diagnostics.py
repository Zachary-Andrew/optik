"""How much of a failed command's output the bench scripts quote.

Two limits cover every case, so a failure reads the same whichever script reports it.
"""

# A check's one-line detail quotes the last DETAIL_CHARS characters of stderr or output.
DETAIL_CHARS = 300

# Reports, raised exceptions and measure.Result keep LOG_TAIL_CHARS, enough for a compiler error or an allocator's message.
# published/common.py scans this tail of measure.Result.err for out-of-memory signatures.
LOG_TAIL_CHARS = 2000
