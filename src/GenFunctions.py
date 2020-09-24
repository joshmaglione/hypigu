#
#   Copyright 2020 Joshua Maglione 
#
#   Distributed under MIT License
#


def UniversalGeneratingFunction(A):
    pass

def LocalIgusaZetaFunction(A):
    pass

def CombinatorialSkeleton(A, Map=False):
    from sage.all import PolynomialRing, QQ, var
    from Globals import __DEFAULT_t as t
    from PosetOps import CharacteristicFunction, PoincarePolynomial
    _, P = CharacteristicFunction(A)
    central = A.is_central()
    if central: 
        prop = filter(lambda X: X != P.top() and X != '', P._elements)
        n = len(prop) + 1
    else: 
        prop = filter(lambda X: X != '', P._elements)
        n = len(prop) 
    P_prop = P.subposet(prop)
    X = PolynomialRing(QQ, 'X', n).gens()
    Y = lambda Z: X[P._elements.index(Z) - 1]
    skele = 0
    for C in P_prop.chains():
        F = [''] + C 
        pi = PoincarePolynomial(A, F)
        Z_in = reduce(lambda x, y: x*Y(y), C, 1)
        C_comp = filter(lambda z: not z in C, P_prop._elements)
        Z_out = reduce(lambda x, y: x*(1 - Y(y)), C_comp, 1)
        skele += pi*Z_in*Z_out
    Skeleton = skele/reduce(lambda x, y: x*(1 - Y(y)), P._elements[1:], 1)
    if Map:
        from sage.all import ZZ
        def user_map(x):
            if x in X:
                s = P._elements[X.index(x) + 1]
            else:
                try:
                    s = P._elements[x + 1]
                except TypeError:
                    raise TypeError("Input must be an integer or variable.")
            return list(map(lambda x: ZZ(int(x)), s.split(' ')))
        return Skeleton, user_map
    return Skeleton