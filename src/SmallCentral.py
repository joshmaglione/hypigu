#
#   Copyright 2020 Joshua Maglione 
#
#   Distributed under MIT License
#

def _small_central(A):
    assert A.is_central()
    P = A.intersection_poset()
    assert P.rank() <= 2
    if P.rank() == 1:
        from Braid import BraidArrangement
        return BraidArrangement(1)
    # Now we assume the rank == 2.
    m = len(P) - 2
    from Globals import __DEFAULT_p, __DEFAULT_t 
    from sage.all import var
    p = var(__DEFAULT_p)
    t = var(__DEFAULT_t)
    return (1 - p**-1)*(1 - (m - 1)*p**-1 + (m - 1)*p**-1*t - p**-2*t) / ((1 - p**-1*t)*(1 - p**-2*t**m))
    