#
#   Copyright 2020 Joshua Maglione
#
#   Distributed under MIT License
#
from sage.all import Partitions, prod, QQ, polygens, PolynomialRing, ZZ, multinomial, factorial, binomial

_TABLE_CUTOFF = 3
_BRAID_TABLE = {
    'Igusa' : [
        lambda p, t: p.parent().one(),
        lambda p, t: (1 - p**-1)/(1 - p**-1*t),
        lambda p, t: (1-p**-1)*(1-2*p**-1+2*p**-1*t-p**-2*t)/((1-p**-1*t)*(1-p**-2*t**3)),
        lambda p, t: p**-6*(1-p**-1)*(p**4*(6-5*p+p**2)-4*p**4*(2-p)*t-p**2*(3-7*p+2*p**2)*t**2+p**2*(2-7*p+3*p**2)*t**3-4*p*(1-2*p)*t**4-(1-5*p+6*p**2)*t**5)/((1-p**-1*t)**2*(1-p**-2*t**3)*(1-p**-3*t**6))
    ],
    'cfHP' : [
        lambda p, t: p.parent().one(),
        lambda p, t: (1 + p)/(1 - t),
        lambda p, t: ((1 + 3*p + 2*p**2) + (2 + 3*p + p**2)*t)/((1 - t)**2),
        lambda p, t: ((1 + 6*p + 11*p**2 + 6*p**3) + (11 + 37*p + 37*p**2 + 11*p**3)*t + (6 + 11*p + 6*p**2 + p**3)*t**2)/((1 - t)**3)
    ],
    'cfHP-red' : [
        lambda t: t.parent().one(),
        lambda t: 2/(1 - t),
        lambda t: (6 + 6*t)/((1 - t)*(1 - t**3)),
        lambda t: (24 + 48*t + 24*t**2 + 48*t**3 + 24*t**4)/((1 - t)*(1 - t**3)*(1 - t**6))
    ]
}

def _Poincare_prod(L):
    """
    Given a partitions L with n parts, construct the product of Poincaré
    polynomials (in Y) of the braid arrangement in ``n-1`` affine space.

    EXAMPLES::

        sage: from hypigu.braid import _Poincare
        sage: _Poincare([1, 1])
        Y + 1
        sage: _Poincare([4, 1, 2])
        2*Y^2 + 3*Y + 1
    """
    A = PolynomialRing(ZZ, 'Y')
    return A.prod(A([1, z]) for z in range(1, len(L)))

def _recursive_crank(p, t, n, known=None, style="Igusa"):
    """
    Constructs the Igusa integral for the braid arrangement
    with little repetition of work.
    """
    # Some functions to use.
    # 1. Counts the number of set partitions of |L| with shape L.
    P = lambda L: multinomial(L) // prod(
        factorial(list(L).count(z)) for z in set(L)
    )
    # 2. Counts the number of edges in a disjoint union of complete graphs
    # determined by L.
    binom_sum = lambda L: sum(binomial(z, 2) for z in L)

    if known is None:
        known = [_BRAID_TABLE[style][k](p, t) for k in range(_TABLE_CUTOFF + 1)]
    k = len(known) + 1
    Zk = 0
    for L in Partitions(k):
        if len(L) > 1:
            if style == "reduced":
                L_factors = [P(L), 1, t**binom_sum(L), factorial(len(L))]
            else:
                L_factors = [
                    P(L), p**(1 - sum(L) + len(L)),
                    t**binom_sum(L), _Poincare_prod(L)(Y=-p**-1)
                ]
            lower_integrals = [known[z - 1] for z in list(L)]
            L_term = prod(L_factors + lower_integrals)
            Zk += L_term
    if style == "reduced":
        Zk /= (1 - t**binomial(k, 2))
    else:
        Zk /= (1 - p**(-k + 1) * t**binomial(k, 2))
    known.append(Zk)
    if k - 1 < n:
        return _recursive_crank(p, t, n, known=known, style=style)
    else:
        return known[n]

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

        sage: Z = BraidArrangementIgusa(2)
        sage: Z
        -(2*t/q - 2/q - t/q^2 + 1)*(1/q - 1)/((t^3/q^2 - 1)*(t/q - 1))
    """
    if n < 0:
        raise ValueError("n must be nonnegative.")
    p, t = polygens(QQ, 'q,t')
    if n <= _TABLE_CUTOFF:
        return _BRAID_TABLE['Igusa'][n](p, t)
    return _recursive_crank(p, t, n)