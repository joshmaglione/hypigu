#
#   Copyright 2020 Joshua Maglione 
#
#   Distributed under MIT License
#

from datetime import datetime as _dt

__PRINT = True
__SANITY = True

def __TIME(): 
    return "[{0}] ".format(_dt.now().strftime("%b %d %H:%M:%S"))