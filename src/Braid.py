#
#   Copyright 2020 Joshua Maglione 
#
#   Distributed under MIT License
#

from sage.all import binomial as _binomial
from sage.all import factorial as _factorial

_TABLE_CUTOFF = 1

def _Igusa_braid_table(p, t, n, style="standard"):
    if n <= 0:
        return 1
    if n == 1:
        if style == "reduced":
            return 2/(1 - t)
        elif style == "skeleton":
            return (1 + p)/(1 - t)
        else:
            return (1 - p**-1)/(1 - p**-1*t)
    if n == 2: 
        if style == "reduced":
            return (6 + 6*t)/((1 - t)*(1 - t**3))
        elif style == "skeleton":
            return ((1 + 3*p + 2*p**2) + (2 + 3*p + p**2)*t)/((1 - t)**2)
        else:
            return (1-p**-1)*(1-2*p**-1+2*p**-1*t-p**-2*t)/((1-p**-1*t)*(1-p**-2*t**3))
    if n == 3: 
        if style == "reduced":
            return (24 + 48*t + 24*t**2 + 48*t**3 + 24*t**4)/((1 - t)*(1 - t**3)*(1 - t**6))
        elif style == "skeleton":
            return ((1 + 6*p + 11*p**2 + 6*p**3) + (11 + 37*p + 37*p**2 + 11*p**3)*t + (6 + 11*p + 6*p**2 + p**3)*t**2)/((1 - t)**3)
        else:
            return p**-6*(1-p**-1)*(p**4*(6-5*p+p**2)-4*p**4*(2-p)*t-p**2*(3-7*p+2*p**2)*t**2+p**2*(2-7*p+3*p**2)*t**3-4*p*(1-2*p)*t**4-(1-5*p+6*p**2)*t**5)/((1-p**-1*t)**2*(1-p**-2*t**3)*(1-p**-3*t**6))
    raise ValueError("Bad n-value.")

# Counts the number of set partitions of [L[0] + ... + L[-1]]. So instead of
# running through all *labeled* partitions, we just run through the partition
# "shape" and count the number of different labelings.
def _P(L):
    from sage.all import Set
    # Compute binom(n, p_1)*binom(n - p_1, p_2)*binom(n - p_1 - p_2, p_3)*...
    def binom(n, P):
        if len(P) == 1:
            return _binomial(n, P[0])
        else:
            return _binomial(n, P[0]) * binom(n - P[0], P[1:])
    n = reduce(lambda x, y: x + y, L)
    count = lambda n: len(filter(lambda x: x == n, L))
    S = list(Set(list(L)))
    d = reduce(lambda x, y: x*y, map(lambda z: _factorial(count(z)), S))
    return binom(n, L) // d

# Counts the number of edges in a subgraph of the complete graph determined by L
def _binom_sum(L):
    binomials = map(lambda z: _binomial(z, 2), L)
    return reduce(lambda x, y: x + y, binomials)

# Constructs the Poincare polynomial (in Y) of the braid arrangement in |L|-1
# affine space.
def _Poincare(L): 
    from sage.all import var
    Y = var('Y')
    factors= map(lambda z: 1 + z*Y, range(1, len(L)))
    return reduce(lambda x, y: x*y, factors, 1)

# Constructs the Igusa integral for the braid arrangement with little repetition
# of work. 
def _recursive_crank(p, t, n, known=[], style="standard"):
    from sage.all import Partitions 
    if known == []:
        known = [_Igusa_braid_table(p, t, k, style=style) 
            for k in range(_TABLE_CUTOFF + 1)]
    k = len(known) + 1
    Zk = 0
    for L in Partitions(k):
        if len(L) > 1:
            if style == "reduced":
                L_factors = [_P(L), 1, t**(_binom_sum(L)), _factorial(len(L))]
            else:
                L_factors = [
                    _P(L), p**(1-reduce(lambda x, y: x + y - 1, L)), 
                    t**(_binom_sum(L)), _Poincare(L)(Y=-p**-1)
                ]
            lower_integrals = map(lambda z: known[z - 1], list(L))
            L_term = reduce(lambda x, y: x*y, L_factors + lower_integrals, 1)
            Zk += L_term
    if style == "reduced":
        Zk = Zk/(1 - t**(_binomial(k, 2)))
    else:
        Zk = Zk/(1 - p**(-k + 1)*t**(_binomial(k, 2)))
    known += [Zk]
    if k - 1 < n:
        return _recursive_crank(p, t, n, known=known, style=style)
    else:
        return known[n]

def BraidArrangementIgusa(n, style="standard"):
    r"""
    Return the rational function associated to the local Igusa zeta function for
    the n-dimensional essential braid arrangement.

    INPUT:

    - ``n`` -- integer; the dimension of the ambient affine space.

    - ``style`` -- string (default: `standard`); the style of computation. 

    OUTPUT: A rational function in at most two variables.

    EXAMPLES:

    This example illustrates what to expect ::

        sage: Z = BraidArrangementIgusa(2)
        sage: Z
        -(2*t/q - 2/q - t/q^2 + 1)*(1/q - 1)/((t^3/q^2 - 1)*(t/q - 1))

    """
    from sage.all import var
    from Globals import __DEFAULT_p, __DEFAULT_t
    if not isinstance(style, str):
        raise TypeError("Expected ``style`` to be a string.")
    if not style.lower() in {"standard", "reduced"}:
        raise ValueError("Unknown style")
    p = var(__DEFAULT_p)
    t = var(__DEFAULT_t)
    if n <= _TABLE_CUTOFF:
        return _Igusa_braid_table(p, t, n, style=style.lower())
    return _recursive_crank(p, t, n, style=style)

