#
#   Copyright 2020 Joshua Maglione
#
#   Distributed under MIT License
#
from functools import reduce

from sage.arith.all import factorial, binomial
from sage.all import Partitions, Set, prod
from sage.all import QQ, polygens, PolynomialRing, ZZ
from sage.arith.misc import multinomial


_TABLE_CUTOFF = 3


def _Igusa_braid_table(p, t, n, style="standard"):
    if n <= 0:
        return p.parent().one()
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
    raise ValueError("bad n-value")


def _P(L):
    """
    Counts the number of set partitions of [L[0] + ... + L[-1]].

    So instead of running through all *labeled* partitions, we just
    run through the partition "shape" and count the number of
    different labelings.

    EXAMPLES::

        sage: from hypigu.Braid import _P
        sage: _P([4,5,6])
        ?
    """
    return multinomial(L) // prod(factorial(L.count(z)) for z in Set(L))


def _binom_sum(L):
    """
    Counts the number of edges in a subgraph of the complete graph
    determined by L.

    EXAMPLES::

        sage: from hypigu.Braid import _binom_sum
        sage: _binom_sum([4,5,6])
        ?
    """
    return sum(binomial(z, 2) for z in L)


def _Poincare(L):
    """
    Construct the Poincaré polynomial (in Y) of the braid arrangement
    in ``|L|-1`` affine space.

    EXAMPLES::

        sage: from hypigu.Braid import _Poincare
        sage: _Poincare(4)
        ?
    """
    A = PolynomialRing(ZZ, 'Y')
    return A.prod(A([1, z]) for z in range(1, len(L)))


def _recursive_crank(p, t, n, known=[], style="standard"):
    """
    Constructs the Igusa integral for the braid arrangement
    with little repetition of work.
    """
    if not known:
        known = [_Igusa_braid_table(p, t, k, style=style)
            for k in range(_TABLE_CUTOFF + 1)]
    k = len(known) + 1
    Zk = 0
    for L in Partitions(k):
        if len(L) > 1:
            if style == "reduced":
                L_factors = [_P(L), 1, t**_binom_sum(L), factorial(len(L))]
            else:
                L_factors = [
                    _P(L), p**(1 - reduce(lambda x, y: x + y - 1, L)),
                    t**_binom_sum(L), _Poincare(L)(Y=-p**-1)
                ]
            lower_integrals = [known[z - 1] for z in list(L)]
            L_term = prod(L_factors + lower_integrals)
            Zk += L_term
    if style == "reduced":
        Zk /= (1 - t**binomial(k, 2))
    else:
        Zk /= (1 - p**(-k + 1) * t**binomial(k, 2))
    known += [Zk]
    if k - 1 < n:
        return _recursive_crank(p, t, n, known=known, style=style)
    else:
        return known[n]


def BraidArrangementIgusa(n):
    r"""
    Return the rational function associated to the local Igusa zeta function for
    the n-dimensional essential braid arrangement.

    INPUT:

    - ``n`` -- integer; the dimension of the ambient affine space.

    OUTPUT: A rational function in at most two variables.

    EXAMPLES:

    This example illustrates what to expect ::

        sage: Z = BraidArrangementIgusa(2)
        sage: Z
        -(2*t/q - 2/q - t/q^2 + 1)*(1/q - 1)/((t^3/q^2 - 1)*(t/q - 1))
    """
    p, t = polygens(QQ, 'q,t')
    if n <= _TABLE_CUTOFF:
        return _Igusa_braid_table(p, t, n, style="standard")
    return _recursive_crank(p, t, n, style="standard")
