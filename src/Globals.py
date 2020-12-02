#
#   Copyright 2020 Joshua Maglione 
#
#   Distributed under MIT License
#

from datetime import datetime as _dt

__DEFAULT_p = 'q'
__DEFAULT_t = 't'
__PRINT = True

def __TIME(): 
    return "[{0}] ".format(_dt.now().strftime("%H:%M:%S"))