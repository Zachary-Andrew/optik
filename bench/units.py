"""Unit conversions shared by the bench scripts.

Memory is measured in binary units (KiB, MiB, GiB, TiB), matching what RSS and cgroup limits report.
Disk and download sizes use decimal units (MB, GB), matching df and curl.
"""

BYTES_PER_KIB = 1024
BYTES_PER_MIB = 1024**2
BYTES_PER_GIB = 1024**3
BYTES_PER_TIB = 1024**4
KIB_PER_MIB = 1024
KIB_PER_GIB = 1024**2
MIB_PER_GIB = 1024

BYTES_PER_MB = 10**6
BYTES_PER_GB = 10**9

BITS_PER_BYTE = 8

NS_PER_SECOND = 1e9
US_PER_SECOND = 1e6
SECONDS_PER_MINUTE = 60

PERCENT = 100
