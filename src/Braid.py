#
#   Copyright 2020 Joshua Maglione 
#
#   Distributed under MIT License
#

from sage.all import binomial as _binomial
from sage.all import factorial as _factorial
from sage.all import Partitions as _partitions
from sage.all import Set as _set
from sage.all import var as _var

_TABLE_CUTOFF = 3

def _Igusa_braid_table(p, t, n):
    if n <= 0:
        return 1
    if n == 1:
        return (1 - p**-1)/(1 - p**-1*t)
    if n == 2: 
        return (1-p**-1)*(1-2*p**-1+2*p**-1*t-p**-2*t)/((1-p**-1*t)*(1-p**-2*t**3))
    if n == 3: 
        return p**-6*(1-p**-1)*(p**4*(6-5*p+p**2)-4*p**4*(2-p)*t-p**2*(3-7*p+2*p**2)*t**2+p**2*(2-7*p+3*p**2)*t**3-4*p*(1-2*p)*t**4-(1-5*p+6*p**2)*t**5)/((1-p**-1*t)**2*(1-p**-2*t**3)*(1-p**-3*t**6))
    raise ValueError("Bad n-value.")

# Counts the number of set partitions of [L[0] + ... + L[-1]]
def _P(L):
    # Compute binom(n, p_1)*binom(n - p_1, p_2)*binom(n - p_1 - p_2, p_3)*...
    def binom(n, P):
        if len(P) == 1:
            return _binomial(n, P[0])
        else:
            return _binomial(n, P[0]) * binom(n - P[0], P[1:])
    n = reduce(lambda x, y: x + y, L)
    count = lambda n: len(filter(lambda x: x == n, L))
    S = list(_set(list(L)))
    d = reduce(lambda x, y: x*y, map(lambda z: _factorial(count(z)), S))
    return binom(n, L) // d

# Counts the number of edges in a subgraph determined by L of the complete graph
def _beta(L):
    binomials = map(lambda z: _binomial(z, 2), L)
    return reduce(lambda x, y: x + y, binomials)

# Constructs the characteristic polynomial of the braid arrangement in |L|-1
# affine space.
def _char(p, L): 
    factors= map(lambda z: p-z, range(1, len(L)))
    return reduce(lambda x, y: x*y, factors, 1)

# Constructs the Igusa integral for the braid arrangement with little repetition
# of work. 
def _recursive_crank(p, t, n, known=[]):
    if known == []:
        known = [_Igusa_braid_table(p, t, k) for k in range(_TABLE_CUTOFF + 1)]
    k = len(known) + 1
    Zk = 0
    for L in _partitions(k):
        if len(L) > 1:
            L_factors = [_P(L), p**(-k + 1), t**(_beta(L)), _char(p, L)]
            lower_integrals = map(lambda z: known[z - 1], list(L))
            L_term = reduce(lambda x, y: x*y, L_factors + lower_integrals, 1)
            Zk += L_term
    Zk = Zk/(1 - p**(-k + 1)*t**(_binomial(k, 2)))
    known += [Zk]
    if k - 1 < n:
        return _recursive_crank(p, t, n, known=known)
    else:
        return known[n]

def BraidArrangement(n):
    p = _var("p")
    t = _var("t")
    if n <= _TABLE_CUTOFF:
        return _Igusa_braid_table(p, t, n)
    return _recursive_crank(p, t, n)

