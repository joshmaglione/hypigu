#
#   Copyright 2020 Joshua Maglione 
#
#   Distributed under MIT License
#

from functools import reduce as _reduce

def _non_Weyl_arrangements(X, n, shift=[0]):
    from sage.all import HyperplaneArrangements, QQ 
    def add_shift(x):
        return list(map(lambda k: [x, k], shift))
    if X == 'I':
        # Not expecting n <= 2 for type I.
        if n == 3:
            return _arrangements_from_roots("A", 2, shift=shift)
        if n == 4: 
            return _arrangements_from_roots("B", 2, shift=shift)
        H = HyperplaneArrangements(QQ, tuple(['x0', 'x1']))
        norms = [tuple([1, 0]), tuple([0, 1]), tuple([1, 1]), tuple([1, -1])]
        cval = 2
        for k in range(n - 4):
            norms.append(tuple([(-1)**(k%2), cval + (k // 2)]))
    else:
        if n == 2:
            return _non_Weyl_arrangements("I", 5, shift=shift)
        print("Type H not yet implemented.")
        return None
    aff_norms = _reduce(lambda x, y: x+y, map(lambda x: add_shift(x), norms),[])
    return H(aff_norms)
    
def _arrangements_from_roots(X, n, shift=[0]):
    from sage.all import HyperplaneArrangements, QQ, RootSystem
    Phi = RootSystem([X, n]).ambient_space()
    d = Phi.dimension()
    H = HyperplaneArrangements(QQ, tuple(['x' + str(i) for i in range(d)]))
    norms = Phi.positive_roots()
    def add_shift(x):
        norm = tuple(map(lambda i: QQ(i), x._vector_().list()))
        return list(map(lambda k: [norm, k], shift))
    aff_norms = _reduce(lambda x, y: x+y, map(lambda x: add_shift(x), norms),[])
    return H(aff_norms)

# Verifies that the Coxeter-theoretic data is expected.
def _Coxeter_check(X, n):
    if n <= 0:
        raise ValueError("Expected a positive rank.")
    if not X in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']:
        raise ValueError(
            "Expected first character to be in {A, B, C, D, E, F, G, H, I}."
        )
    if X in ['E', 'F', 'G', 'H']:
        if X == 'E' and not n in {6, 7, 8}:
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
        if X == 'H' and not n in {2, 3, 4}:
            raise ValueError(
                "Expected rank to be either 2, 3, or 4."
            )
        if X == 'I' and n <= 2:
            raise ValueError(
                "Expected rank above 2."
            )
    return True

# Parses the Coxeter type input and runs the check for good input. 
def _parse_Coxeter_input(name, n):
    from sage.all import ZZ
    if not isinstance(name, str):
        try:
            n = ZZ(name)
            name = "A"
        except:
            raise TypeError("Expected ``name`` to be a string.")
    if len(name) >= 2 and n == 0:
        try:
            n = int(name[1:])
        except ValueError:
            print("Expected second character to be an integer.")
    X = name[0].upper()
    _ = _Coxeter_check(X, n)
    return (X, n)

# Builds the direct sum arrangement of arrangements A and B.
def _direct_sum(A, B):
    from sage.all import HyperplaneArrangements as HA
    K = A.base_ring()
    d = A.dimension() + B.dimension()
    A_H = map(lambda H: H.coefficients(), A.hyperplanes())
    B_H = map(lambda H: H.coefficients(), B.hyperplanes())
    HH = HA(K, tuple(['x' + str(k) for k in range(d)]))
    A_emb = map(lambda L: L + [0]*(B.dimension()), A_H)
    B_emb = map(lambda L: L[0:1] + [0]*(A.dimension()) + L[1:], B_H)
    return HH(list(A_emb) + list(B_emb))

def CoxeterArrangement(name, n=0):
    r"""
    Return the Coxeter arrangement of the prescribed type.

    INPUT:

    - ``name`` -- string; the Coxeter group name. The string must include a 
      letter from {A, B, C, D}, and it can also include an integer.

    - ``n`` -- integer (default: `0`); the rank of the Coxeter group. This is 
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
    X, n = _parse_Coxeter_input(name, n)
    if X in {'I', 'H'}:
        return _non_Weyl_arrangements(X, n)
    if [X, n] == ['D', 1]:
        X, n = ('A', 1)
    return _arrangements_from_roots(X, n)


def ShiArrangement(name, n=0):
    r"""
    Return the Shi arrangement associated to a Coxeter arrangement of the
    prescribed type.

    INPUT:

    - ``name`` -- string; the Coxeter group name. The string must include a 
      letter from {A, B, C, D}, and it can also include an integer.

    - ``n`` -- integer (default: `0`); the rank of the Coxeter group. This is 
      not required if ``name`` includes the rank. 

    OUTPUT: the Shi arrangement given as a hyperplane arrangement.

    EXAMPLES:

    This example illustrates how to build a Shi arrangement ::

        sage: A = ShiArrangement("D", 4)
        sage: A
        Arrangement of 16 hyperplanes of dimension 4 and rank 4

    We can also combine the rank into the string as follows ::

        sage: A = ShiArrangement("D4")
        sage: A
        Arrangement of 16 hyperplanes of dimension 4 and rank 4

    """
    X, n = _parse_Coxeter_input(name, n)
    if X in {'I', 'H'}:
        return _non_Weyl_arrangements(X, n, shift=[0, 1])
    return _arrangements_from_roots(X, n, shift=[0, 1])


def LinialArrangement(name, n=0):
    r"""
    Return the Linial arrangement of the prescribed type.

    INPUT:

    - ``name`` -- string; the Coxeter group name. The string must include a 
      letter from {A, B, C, D}, and it can also include an integer.

    - ``n`` -- integer (default: `0`); the rank of the Coxeter group. This is 
      not required if ``name`` includes the rank. 

    OUTPUT: the Linial arrangement given as a hyperplane arrangement.

    EXAMPLES:

    This example illustrates how to build a Linial arrangement ::

        sage: A = LinialArrangement(3)
        sage: A
        Arrangement of 16 hyperplanes of dimension 4 and rank 4

    """

    X, n = _parse_Coxeter_input(name, n)
    if X in {'I', 'H'}:
        return _non_Weyl_arrangements(X, n, shift=[1])
    return _arrangements_from_roots(X, n, shift=[1])


def CatalanArrangement(name, n=0):
    r"""
    Return the Catalan arrangement of the prescribed type.

    INPUT:

    - ``name`` -- string; the Coxeter group name. The string must include a 
      letter from {A, B, C, D}, and it can also include an integer.

    - ``n`` -- integer (default: `0`); the rank of the Coxeter group. This is 
      not required if ``name`` includes the rank. 

    OUTPUT: the Catalan arrangement given as a hyperplane arrangement.

    EXAMPLES:

    This example illustrates how to build a Catalan arrangement ::

        sage: A = CatalanArrangement("D", 4)
        sage: A
        Arrangement of 36 hyperplanes of dimension 4 and rank 4

    """

    X, n = _parse_Coxeter_input(name, n)
    if X in {'I', 'H'}:
        return _non_Weyl_arrangements(X, n, shift=[-1, 0, 1])
    return _arrangements_from_roots(X, n, shift=[-1, 0, 1])

def DirectSum(*args):
    if len(args) == 1:
        HPAs = args[0]
    else:
        HPAs = args
    if len(HPAs) == 1: 
        return HPAs[0]
    else:
        D = _direct_sum(HPAs[0], HPAs[1])
        return _reduce(lambda A, B: _direct_sum(A, B), HPAs[2:], D)


def PolynomialToArrangement(f):
    from .GenFunctions import _parse_poly
    A, M = _parse_poly(f)
    return A
    