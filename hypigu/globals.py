#
#   Copyright 2020--2024 Joshua Maglione
#
#   Distributed under MIT License
#
from os import cpu_count 

verbose = False
ncpus = max(1, cpu_count())

def my_time():
    from datetime import datetime as dt
    return f"[{dt.now().strftime("%b %d %H:%M:%S")}]"

def my_print(verbose:bool, s:str, level:int=0):
    if verbose:
        print(f"{my_time()}{' ' + '\t'*level}{s}")
    return None
