#
#   Copyright 2020 Joshua Maglione 
#
#   Distributed under MIT License
#

# The complete solutions for small central arrangements of rank <= 2.
def small_central(A):
    assert A.is_central()
    assert A.rank() <= 2
    if A.rank() == 1:
        from Braid import BraidArrangement
        return BraidArrangement(1)
    # Now we assume the rank == 2.
    m = len(A)
    from Globals import __DEFAULT_p, __DEFAULT_t 
    from sage.all import var
    p = var(__DEFAULT_p)
    t = var(__DEFAULT_t)
    return (1 - p**-1)*(1 - (m - 1)*p**-1 + (m - 1)*p**-1*t - p**-2*t) / ((1 - p**-1*t)*(1 - p**-2*t**m))

# The direct version of the combinatorial skeleton.
def _comb_skele_direct(A):
    from sage.all import PolynomialRing, QQ, var
    from Globals import __DEFAULT_t, __DEFAULT_p
    from PosetOps import CharacteristicFunction, PoincarePolynomial, _proper_part
    _, P = CharacteristicFunction(A)
    P_prop = _proper_part(P)
    n = len(P_prop)
    t = var(__DEFAULT_t)
    skele = 0
    for C in P_prop.chains():
        F = [''] + C 
        pi = PoincarePolynomial(A, F)
        Z_in = t**(len(C))
        C_comp = filter(lambda z: not z in C, P_prop._elements)
        Z_out = (1 - t)**(len(C_comp))
        skele += pi*Z_in*Z_out
    return skele/((1 - t)**n)

# The direct version of the local Igusa zeta function computation.
def _local_Igusa_direct(A):
    from sage.all import PolynomialRing, QQ, var
    from Globals import __DEFAULT_t, __DEFAULT_p
    from PosetOps import CharacteristicFunction, PoincarePolynomial, _proper_part
    p = var(__DEFAULT_p)
    t = var(__DEFAULT_t)
    _, P = CharacteristicFunction(A)
    P_prop = _proper_part(P)
    hypers = list(filter(lambda x: P.covers('', x), P))
    nHypers = lambda Z: len(list(filter(lambda H: P.le(H, Z), hypers)))
    Y = lambda Z: p**(-len(Z.split(' ')))*t**nHypers(Z)
    skele = 0
    for C in P_prop.chains():
        F = [''] + C 
        pi = PoincarePolynomial(A, F)
        Z_in = reduce(lambda x, y: x*Y(y), C, 1)
        C_comp = filter(lambda z: not z in C, P_prop._elements)
        Z_out = reduce(lambda x, y: x*(1 - Y(y)), C_comp, 1)
        skele += pi.subs({p: -p**-1})*Z_in*Z_out
    return skele/reduce(lambda x, y: x*(1 - Y(y)), P._elements[1:], 1)

# The direct version of the universal generating function computation.
def _universal_gen_func(A, MAP):
    from sage.all import PolynomialRing, QQ, var, ZZ
    from Globals import __DEFAULT_t, __DEFAULT_p
    from PosetOps import CharacteristicFunction, PoincarePolynomial, _proper_part
    p = var(__DEFAULT_p)
    Yvar = var('Y')
    _, P = CharacteristicFunction(A)
    P_prop = _proper_part(P)
    n = len(P_prop)
    X = PolynomialRing(QQ, 'X', n).gens()
    Y = lambda Z: X[P._elements.index(Z) - 1]
    skele = 0
    for C in P_prop.chains():
        F = [''] + C 
        pi = PoincarePolynomial(A, F)
        Z_in = reduce(lambda x, y: x*Y(y), C, 1)
        C_comp = filter(lambda z: not z in C, P_prop._elements)
        Z_out = reduce(lambda x, y: x*(1 - Y(y)), C_comp, 1)
        skele += pi.subs({p: Yvar})*Z_in*Z_out
    Skeleton = skele/reduce(lambda x, y: x*(1 - Y(y)), P._elements[1:], 1)
    if MAP:
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


def CombinatorialSkeleton(A, method="direct"):
    if method == "direct":
        return _comb_skele_direct(A)
    else:
        return None

def LocalIgusaZetaFunction(A, method="direct"):
    if method == "direct":
        return _local_Igusa_direct(A)
    else:
        return None

def UniversalGeneratingFunction(A, Map=False, method="direct"):
    if method == "direct":
        return _universal_gen_func(A, Map)
    else:
        return None