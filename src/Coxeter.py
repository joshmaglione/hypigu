#
#   Copyright 2020 Joshua Maglione 
#
#   Distributed under MIT License
#

from sage.all import hyperplane_arrangements as _defaults
from sage.all import HyperplaneArrangements as _HA
from sage.all import QQ as _QQ

def _build_hyperplanes(X, n):
    if X == "A":
        return _defaults.braid(n + 1)
    if X in ["B", "C", "D"]:
        H = _HA(_QQ, tuple(['x' + str(i) for i in range(n)]))
        varbs = H.gens()
        H_ij = [varbs[i] - varbs[j] 
            for i in range(n) for j in range(i + 1, n)]
        H_ij_m = [varbs[i] + varbs[j] 
            for i in range(n) for j in range(i + 1, n)]
        H_k = [varbs[k] for k in range(n)]
        if X == "D": 
            return H(H_ij + H_ij_m)
        else: 
            return H(H_ij + H_ij_m + H_k)


def CoxeterArrangement(name, n=0):
    r"""
    Return the Coxeter arrangement of the prescribed type.

    INPUT:

    - ``name`` -- string; the Coxeter group name. The string must include a 
      letter from {A, B, C, D}, and it can also include an integer.

    - ``y`` -- integer (default: `0`); the rank of the Coxeter group. This is 
      not required if ``name`` includes the rank. 

    OUTPUT: the Coxeter arrangement given as a hyperplane arrangement.

    EXAMPLES:

    This example illustrates how to build a Coxeter arrangement ::

        sage: A = CoxeterArrangement("B", 4)
        sage: A
        Arrangement of 16 hyperplanes of dimension 4 and rank 4

    We can also combine the rank into the string as follows ::

        sage: A = IA.CoxeterArrangement("B4")
        sage: A
        Arrangement of 16 hyperplanes of dimension 4 and rank 4

    """
    if not isinstance(name, str):
        raise TypeError("Expected ``name`` to be a string.")
    if len(name) >= 2 and n == 0:
        try:
            n = int(name[1])
        except ValueError:
            print("Expected second character to be an integer.")
    if n == 0:
        raise ValueError("Expected a positive rank.")
    if not name[0].upper() in ['A', 'B', 'C', 'D']:
        raise ValueError("Expected first character to be in {A, B, C, D}.")
    else:
        X = name[0].upper()
    return _build_hyperplanes(X, n)