#
#   Copyright 2020--2024 Joshua Maglione
#
#   Distributed under MIT License
#

from datetime import datetime as dt
from os import cpu_count 

_PRINT = False
_SANITY_CHECK = False
_NCPUS = max(1, cpu_count() - 1)


def _TIME():
    return f"[{dt.now().strftime("%b %d %H:%M:%S")}] "
