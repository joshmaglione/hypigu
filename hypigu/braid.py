#
#   Copyright 2020--2024 Joshua Maglione
#
#   Distributed under MIT License
#
from sage.combinat.partition import Partitions
from sage.misc.misc_c import prod
from sage.rings.rational_field import Q as QQ
from sage.rings.polynomial.polynomial_ring import polygens
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.integer_ring import Z as ZZ
from sage.arith.misc import multinomial
from sage.functions.other import factorial
from sage.functions.other import binomial

_TABLE_CUTOFF = 5
_BRAID_TABLE = [
    lambda p, t: p.parent().one(),
    lambda p, t: (1 - p)/(1 - p*t),
    lambda p, t: (1-p)*(1-2*p+2*p*t-p**2*t)/((1-p*t)*(1-p**2*t**3)),
    lambda p, t: (p - 1) * (p*t - 1)**-2 * (p*t**2 - 1)**-1 * (p**2*t**3 - 1)**-1 * (p**2*t**4 + p*t**2 + 1)**-1 * (p**6*t**5 - 5*p**5*t**5 + 4*p**5*t**4 + 6*p**4*t**5 - 8*p**4*t**4 - 2*p**4*t**3 + 3*p**4*t**2 + 7*p**3*t**3 - 7*p**3*t**2 - 3*p**2*t**3 + 2*p**2*t**2 + 8*p**2*t - 6*p**2 - 4*p*t + 5*p - 1),
    lambda p, t: (p - 1) * (p*t - 1)**-2 * (p*t**2 - 1)**-1 * (p**2*t**3 - 1)**-1 * (p**2*t**4 + p*t**2 + 1)**-1 * (p**2*t**5 - 1)**-1 * (p**2*t**5 + 1)**-1 * (p**10*t**11 - 9*p**9*t**11 + 8*p**9*t**10 + 26*p**8*t**11 - 42*p**8*t**10 - 24*p**7*t**11 + 6*p**8*t**9 + 58*p**7*t**10 + 9*p**8*t**8 - 9*p**7*t**9 - 12*p**6*t**10 - 41*p**7*t**8 - 9*p**6*t**9 + 12*p**7*t**7 + 54*p**6*t**8 + 6*p**5*t**9 - 18*p**6*t**7 - 16*p**5*t**8 + 4*p**7*t**5 - 6*p**6*t**6 - 18*p**5*t**7 - 21*p**6*t**5 + 29*p**5*t**6 + 12*p**4*t**7 + 12*p**6*t**4 + 29*p**5*t**5 - 21*p**4*t**6 - 18*p**5*t**4 - 6*p**4*t**5 + 4*p**3*t**6 - 16*p**5*t**3 - 18*p**4*t**4 + 6*p**5*t**2 + 54*p**4*t**3 + 12*p**3*t**4 - 9*p**4*t**2 - 41*p**3*t**3 - 12*p**4*t - 9*p**3*t**2 + 9*p**2*t**3 + 58*p**3*t + 6*p**2*t**2 - 24*p**3 - 42*p**2*t + 26*p**2 + 8*p*t - 9*p + 1),
    lambda p, t: (p - 1) * (p*t - 1)**-3 * (p*t**2 - 1)**-1 * (p*t**3 - 1)**-1 * (p**2*t**3 - 1)**-2 * (p**2*t**4 + p*t**2 + 1)**-1 * (p**2*t**5 - 1)**-1 * (p**2*t**5 + 1)**-1 * (p**4*t**12 + p**3*t**9 + p**2*t**6 + p*t**3 + 1)**-1 * (p**18*t**25 - 14*p**17*t**25 + 12*p**17*t**24 + 71*p**16*t**25 - 108*p**16*t**24 - 154*p**15*t**25 + 18*p**16*t**23 + 312*p**15*t**24 + 120*p**14*t**25 + 18*p**16*t**22 - 57*p**15*t**23 - 288*p**14*t**24 - 148*p**15*t**22 - 72*p**14*t**23 + 36*p**15*t**21 + 422*p**14*t**22 + 273*p**13*t**23 - 84*p**14*t**21 - 428*p**13*t**22 - 90*p**12*t**23 + 14*p**15*t**19 - 36*p**14*t**20 - 204*p**13*t**21 + 64*p**12*t**22 - 115*p**14*t**19 + 294*p**13*t**20 + 516*p**12*t**21 + 48*p**14*t**18 + 312*p**13*t**19 - 456*p**12*t**20 - 120*p**11*t**21 - 195*p**13*t**18 - 453*p**12*t**19 - 6*p**11*t**20 - 48*p**13*t**17 + 150*p**12*t**18 + 418*p**11*t**19 + 60*p**10*t**20 + 5*p**14*t**15 + 12*p**13*t**16 + 330*p**12*t**17 + 255*p**11*t**18 - 128*p**10*t**19 - 46*p**13*t**15 + 28*p**12*t**16 - 705*p**11*t**17 - 378*p**10*t**18 + 30*p**13*t**14 + 43*p**12*t**15 - 338*p**11*t**16 + 570*p**10*t**17 + 120*p**9*t**18 - 132*p**12*t**14 + 388*p**11*t**15 + 627*p**10*t**16 - 147*p**9*t**17 - 30*p**12*t**13 + 54*p**11*t**14 - 732*p**10*t**15 - 360*p**9*t**16 + 30*p**12*t**12 + 231*p**11*t**13 - 36*p**10*t**14 + 114*p**9*t**15 + 89*p**8*t**16 - 128*p**11*t**12 - 326*p**10*t**13 + 684*p**9*t**14 + 60*p**8*t**15 - 10*p**7*t**16 - 60*p**11*t**11 + 46*p**10*t**12 - 71*p**9*t**13 - 444*p**8*t**14 + 444*p**10*t**11 + 71*p**9*t**12 - 46*p**8*t**13 + 60*p**7*t**14 + 10*p**11*t**9 - 60*p**10*t**10 - 684*p**9*t**11 + 326*p**8*t**12 + 128*p**7*t**13 - 89*p**10*t**9 - 114*p**9*t**10 + 36*p**8*t**11 - 231*p**7*t**12 - 30*p**6*t**13 + 360*p**9*t**9 + 732*p**8*t**10 - 54*p**7*t**11 + 30*p**6*t**12 + 147*p**9*t**8 - 627*p**8*t**9 - 388*p**7*t**10 + 132*p**6*t**11 - 120*p**9*t**7 - 570*p**8*t**8 + 338*p**7*t**9 - 43*p**6*t**10 - 30*p**5*t**11 + 378*p**8*t**7 + 705*p**7*t**8 - 28*p**6*t**9 + 46*p**5*t**10 + 128*p**8*t**6 - 255*p**7*t**7 - 330*p**6*t**8 - 12*p**5*t**9 - 5*p**4*t**10 - 60*p**8*t**5 - 418*p**7*t**6 - 150*p**6*t**7 + 48*p**5*t**8 + 6*p**7*t**5 + 453*p**6*t**6 + 195*p**5*t**7 + 120*p**7*t**4 + 456*p**6*t**5 - 312*p**5*t**6 - 48*p**4*t**7 - 516*p**6*t**4 - 294*p**5*t**5 + 115*p**4*t**6 - 64*p**6*t**3 + 204*p**5*t**4 + 36*p**4*t**5 - 14*p**3*t**6 + 90*p**6*t**2 + 428*p**5*t**3 + 84*p**4*t**4 - 273*p**5*t**2 - 422*p**4*t**3 - 36*p**3*t**4 + 72*p**4*t**2 + 148*p**3*t**3 + 288*p**4*t + 57*p**3*t**2 - 18*p**2*t**3 - 120*p**4 - 312*p**3*t - 18*p**2*t**2 + 154*p**3 + 108*p**2*t - 71*p**2 - 12*p*t + 14*p - 1)
]

def _recursive_crank(p, t, n, known=None):
    """
    Constructs the Igusa integral for the braid arrangement
    with little repetition of work.
    """
    # 1. Counts the number of set partitions of |L| with shape L.
    P = lambda L: multinomial(L) // prod(
        factorial(list(L).count(z)) for z in set(L)
    )
    # 2. Counts the number of edges in a disjoint union of complete graphs
    #    determined by L.
    binom_sum = lambda L: sum(binomial(z, 2) for z in L)
    # 3. Constructs the product of Poincaré polynomials for the braid
    #    arrangement
    def Poincare_prod(L): 
        A = PolynomialRing(ZZ, 'Y')
        return A.prod(A([1, z]) for z in range(1, len(L)))(Y=-p)

    if known is None:
        known = [_BRAID_TABLE[k](p, t) for k in range(_TABLE_CUTOFF + 1)]
    k = len(known) + 1
    Zk = 0
    for L in (L for L in Partitions(k) if len(L) > 1):
        L_factors = [
            P(L), p**(sum(L) - len(L)), t**binom_sum(L), Poincare_prod(L)
        ]
        lower_integrals = [known[z - 1] for z in list(L)]
        L_term = prod(L_factors + lower_integrals)
        Zk += L_term
    Zk /= (1 - p**(k - 1) * t**binomial(k, 2))
    known.append(Zk)
    if k - 1 < n:
        return _recursive_crank(p, t, n, known=known)
    else:
        return known[n]

# Based on Lemma 5.14 of Maglione--Voll
def BraidArrangementIgusa(n:int):
    r"""
    Return the rational function associated to the local Igusa zeta function for
    the n-dimensional essential braid arrangement. This function is much faster
    than the all-purpose IgusaZetaFunction since it takes advantage of the
    well-known combinatorics of the braid arrangement. 

    INPUT:

    - ``n`` -- integer; the dimension of the ambient affine space.

    OUTPUT: A rational function in variables q and t.

    EXAMPLES:

    This example illustrates what to expect ::

        sage: Z = hi.BraidArrangementIgusa(2)
        sage: Z
        -(2*t/q - 2/q - t/q^2 + 1)*(1/q - 1)/((t^3/q^2 - 1)*(t/q - 1))
    """
    if n < 0:
        raise ValueError("n must be nonnegative.")
    y, t = polygens(QQ, 'y,t')
    if n <= _TABLE_CUTOFF:
        return _BRAID_TABLE[n](y, t)
    return _recursive_crank(y, t, n)