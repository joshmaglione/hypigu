#
#   Copyright 2020--2024 Joshua Maglione
#
#   Distributed under MIT License
#

from functools import reduce
from sage.geometry.hyperplane_arrangement.arrangement import HyperplaneArrangements
from sage.rings.rational_field import Q as QQ
from sage.combinat.root_system.coxeter_group import CoxeterGroup
from sage.combinat.root_system.root_system import RootSystem
from sage.combinat.subset import Subsets
from sage.matrix.constructor import Matrix
from sage.symbolic.ring import SR

# Given a polynomial, return a hyperplane arrangement equivalent to the linear
# factors of f.
def _parse_poly(f):
    if isinstance(f, str):
        f = SR(f)
    if f.base_ring() == SR:
        L = f.factor_list()
        K = QQ
    else:
        L = list(f.factor())
        K = f.base_ring()

    L = filter(lambda T: not T[0] in K, L)  # Remove constant factors
    F, M = list(zip(*L))

    # Verify that each polynomial factor is linear
    is_lin = lambda g: all(g.degree(x) <= 1 for x in g.variables())
    if not all(map(is_lin, F)):
        raise ValueError("Expected product of linear factors.")

    varbs = f.variables()
    varbs_str = tuple(str(x) for x in varbs)
    HH = HyperplaneArrangements(K, varbs_str)

    def poly_vec(g):
        c = K(g.subs({x: 0 for x in g.variables()}))
        return tuple([c] + [K(g.coefficient(x)) for x in varbs])

    F_vec = tuple(map(poly_vec, F))
    A = HH(Matrix(K, F_vec))

    # This scrambles the hyperplanes, so we need to scramble M in the same way.
    A_vec = tuple(map(lambda H: tuple(H.coefficients()), A.hyperplanes()))
    perm = (F_vec.index(v) for v in A_vec)
    M_new = tuple([M[i] for i in perm])

    return A, M_new

def _non_Weyl_arrangements(X, n, shift=[0]):
    if X == 'I':
        # Not expecting n <= 2 for type I.
        if n == 3:
            return _arrangements_from_roots("A", 2, shift=shift)
        if n == 4:
            return _arrangements_from_roots("B", 2, shift=shift)
        H = HyperplaneArrangements(QQ, ('x0', 'x1'))
        norms = [(1, 0), (0, 1), (1, 1), (1, -1)]
        norms += [((-1)**(k % 2), 2 + (k // 2)) for k in range(n - 4)]
    else:
        if n == 2:
            return _non_Weyl_arrangements("I", 5, shift=shift)
        W = CoxeterGroup([X, n])
        H = HyperplaneArrangements(
            W.base_ring(), tuple(['x' + str(k) for k in range(n)])
        )
        norms = W.positive_roots()
    def add_shift(x):
        return [[k] + list(x) for k in shift]
    aff_norms = reduce(lambda x, y: x + add_shift(y), norms, [])
    return H(aff_norms)


def _arrangements_from_roots(X, n, shift=[0]):
    Phi = RootSystem([X, n]).ambient_space()
    d = Phi.dimension()
    H = HyperplaneArrangements(QQ, tuple(['x' + str(i) for i in range(d)]))
    norms = Phi.positive_roots()

    def add_shift(x):
        norm = tuple(map(lambda i: QQ(i), x._vector_().list()))
        return [[norm, k] for k in shift]
    aff_norms = reduce(lambda x, y: x + add_shift(y), norms, [])
    return H(aff_norms)


# Verifies that the Coxeter-theoretic data is expected.
def _Coxeter_check(X, n):
    if n <= 0:
        raise ValueError("Expected a positive rank.")
    if X not in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']:
        raise ValueError(
            "Expected first character to be in {A, B, C, D, E, F, G, H, I}."
        )
    if X in ['E', 'F', 'G', 'H']:
        if X == 'E' and n not in {6, 7, 8}:
            raise ValueError(
                "Expected rank to be either 6, 7, or 8."
            )
        if X == 'F' and n != 4:
            raise ValueError(
                "Expected rank to be 4."
            )
        if X == 'G' and n != 2:
            raise ValueError(
                "Expected rank to be 2."
            )
        if X == 'H' and n not in {2, 3, 4}:
            raise ValueError(
                "Expected rank to be either 2, 3, or 4."
            )
        if X == 'I' and n <= 2:
            raise ValueError(
                "Expected rank above 2."
            )
    return True


# Parses the Coxeter type input and runs the check for good input.
def _parse_Coxeter_input(name):
    if isinstance(name, str):
        factors = name.split(' ')
    else:
        factors = list(name)
        if not isinstance(factors[0], str):
            raise TypeError("Expected ``name`` to be an interable container of strings.")

    def convert(s):
        if len(s) <= 1:
            raise ValueError(f"Cannot parse input: {s}")
        try:
            n = int(s[1:])
        except ValueError:
            raise ValueError(f"Cannot parse as integer: {s[1:]}")
        return (s[0].upper(), n)
    Cox_facts = [convert(s) for s in factors]
    assert all(map(lambda X: _Coxeter_check(X[0], X[1]), Cox_facts))
    return Cox_facts


# Builds the direct sum arrangement of arrangements A and B.
def _direct_sum(A, B):
    K = A.base_ring()
    d = A.dimension() + B.dimension()
    A_H = map(lambda H: H.coefficients(), A.hyperplanes())
    B_H = map(lambda H: H.coefficients(), B.hyperplanes())
    HH = HyperplaneArrangements(K, tuple(['x' + str(k) for k in range(d)]))
    A_emb = map(lambda L: L + [0]*(B.dimension()), A_H)
    B_emb = map(lambda L: L[0:1] + [0]*(A.dimension()) + L[1:], B_H)
    return HH(list(A_emb) + list(B_emb))


def _basic_wrapper(name, s):
    factors = _parse_Coxeter_input(name)
    def irr_facts(X):
        if X[0] in {'I', 'H'}:
            return _non_Weyl_arrangements(X[0], X[1], shift=s)
        if X == ['D', 1]:
            X = ['A', 1]
        return _arrangements_from_roots(X[0], X[1], shift=s)
    irr_arr = [irr_facts(X) for X in factors]
    return DirectSum(irr_arr)


def _res_arr(n):
    Subs = Subsets(range(1, n + 1))
    def S_vec(S):
        v = [QQ(0)]*(n+1)
        for s in S:
            v[s] = QQ(1)
        return v
    M = [S_vec(S) for S in Subs if len(S) > 0]
    HH = HyperplaneArrangements(QQ, tuple(['x' + str(k) for k in range(n)]))
    H = HH(Matrix(QQ, M))
    return H


def DirectSum(*args):
    r"""
    Return the direct sum of hyperplane arrangements.

    INPUT:

    - An iterable container of hyperplane arrangements.

    OUTPUT: the direct sum arrangement given as a hyperplane arrangement.

    EXAMPLE:

        sage: A = hi.CoxeterArrangement("A1")
        sage: Bool3 = hi.DirectSum(A, A, A)
        sage: Bool3
        Arrangement <x4 - x5 | x2 - x3 | x0 - x1>
        sage: Bool64 = hi.DirectSum([A]*64)
        sage: Bool64
        Arrangement of 64 hyperplanes of dimension 128 and rank 64
    """
    if len(args) == 1:
        HPAs = args[0]
    else:
        HPAs = args
    if len(HPAs) == 1:
        return HPAs[0]
    else:
        D = _direct_sum(HPAs[0], HPAs[1])
        return reduce(lambda A, B: _direct_sum(A, B), HPAs[2:], D)


def CoxeterArrangement(name):
    r"""
    Return the Coxeter arrangement of the prescribed Coxeter type.

    INPUT:

    - ``name`` -- an iterable container of strings; the Coxeter group
        name. The string must include a letter from {A, ..., H} and a
        nonnegative integer.

    OUTPUT: the Coxeter arrangement given as a hyperplane arrangement.

    EXAMPLES:

        sage: hi.CoxeterArrangement("B4")
        Arrangement of 16 hyperplanes of dimension 4 and rank 4

        sage: hi.CoxeterArrangement("A1 A2")
        Arrangement <x3 - x4 | x2 - x3 | x2 - x4 | x0 - x1>

        sage: hi.CoxeterArrangement(["D4", "E6"])
        Arrangement of 48 hyperplanes of dimension 12 and rank 10
    """
    return _basic_wrapper(name, [0])


def ShiArrangement(name):
    r"""
    Return the Shi arrangement of the Coxeter prescribed type.

    INPUT:

    - ``name`` -- an iterable container of strings; the Coxeter group
        name. The string must include a letter from {A, ..., H} and a
        nonnegative integer.

    OUTPUT: the Shi arrangement given as a hyperplane arrangement.

    EXAMPLES:

        sage: hi.ShiArrangement("A2")
        Arrangement of 6 hyperplanes of dimension 3 and rank 2

        sage: hi.ShiArrangement("A1 A1")
        Arrangement <x2 - x3 | x2 - x3 + 1 | x0 - x1 | x0 - x1 + 1>

        sage: hi.ShiArrangement(["H3", "B3"])
        Arrangement of 48 hyperplanes of dimension 6 and rank 6
    """
    return _basic_wrapper(name, [0, 1])


def LinialArrangement(name):
    r"""
    Return the Linial arrangement of the prescribed Coxeter type.

    INPUT:

    - ``name`` -- an iterable container of strings; the Coxeter group
        name. The string must include a letter from {A, ..., H} and a
        nonnegative integer.

    OUTPUT: the Linial arrangement given as a hyperplane arrangement.

    EXAMPLES:

        sage: hi.LinialArrangement("D4")
        Arrangement of 12 hyperplanes of dimension 4 and rank 4

        sage: hi.LinialArrangement("A1 A1")
        Arrangement <x2 - x3 + 1 | x0 - x1 + 1>

        sage: hi.LinialArrangement(["I4", "F4"])
        Arrangement of 28 hyperplanes of dimension 6 and rank 6
    """
    return _basic_wrapper(name, [1])


def CatalanArrangement(name):
    r"""
    Return the Catalan arrangement of the prescribed Coxeter type.

    INPUT:

    - ``name`` -- an iterable container of strings; the Coxeter group
        name. The string must include a letter from {A, ..., H} and a
        nonnegative integer.

    OUTPUT: the Catalan arrangement given as a hyperplane arrangement.

    EXAMPLES:

        sage: hi.CatalanArrangement("B4")
        Arrangement of 48 hyperplanes of dimension 4 and rank 4

        sage: hi.CatalanArrangement("A1 A1")
        Arrangement of 6 hyperplanes of dimension 4 and rank 2

        sage: hi.CatalanArrangement(["D4", "A3"])
        Arrangement of 54 hyperplanes of dimension 8 and rank 7
    """
    return _basic_wrapper(name, [-1, 0, 1])


def ResonanceArrangement(n):
    r"""
    Return the rank-n resonance arrangement.

    INPUT:

    - ``n`` -- a positive integer.

    OUTPUT: the resonance arrangement given as a hyperplane arrangement.

    EXAMPLES:

        sage: hi.ResonanceArrangement(3)
        Arrangement of 7 hyperplanes of dimension 3 and rank 3

        sage: hi.ResonanceArrangement(2)
        Arrangement <x1 | x0 | x0 + x1>

        sage: hi.ResonanceArrangement(10)
        Arrangement of 1023 hyperplanes of dimension 10 and rank 10
    """
    return _res_arr(n)


def PolynomialToArrangement(f):
    r"""
    Return the associated hyperplane arrangement of f, the given product of
    linear polynomials.

    INPUT:

    - ``f`` -- a polynomial or symbolic expression.

    OUTPUT: the hyperplane arrangement associated with the linear factors of f.

    EXAMPLES:

        sage: f = SR('X^3*Y^2*Z')
        sage: f
        X^3*Y^2*Z
        sage: hi.PolynomialToArrangement(f)
        Arrangement <Z | Y | X>

        sage: f = 'X*Y*Z*W*(X - Y)*(X - Z)*(Y - W)'
        sage: hi.PolynomialToArrangement(f)
        Arrangement of 7 hyperplanes of dimension 4 and rank 4
    """
    A, _ = _parse_poly(f)
    return A