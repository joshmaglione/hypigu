#
#   Copyright 2020 Joshua Maglione
#
#   Distributed under MIT License
#

from datetime import datetime as _dt
from os import cpu_count as _cpu_count

__PRINT = False
__SANITY = False
__NCPUS = max(1, _cpu_count()-1)

def __TIME(): 
    return "[{0}] ".format(_dt.now().strftime("%b %d %H:%M:%S"))