#
#   Copyright 2020 Joshua Maglione 
#
#   Distributed under MIT License
#

from Database import internal_database as _data

# The complete solutions for small central arrangements of rank <= 2.
def _small_central(A, style):
    from Globals import __DEFAULT_p, __DEFAULT_t 
    from sage.all import var
    assert A.is_central()
    assert A.rank() <= 2
    p = var(__DEFAULT_p)
    t = var(__DEFAULT_t)
    if A.rank() == 1:
        if style == 'Igusa':
            return (1 - p**-1)/(1 - p**-1*t)
        if style == 'skele':
            return (1 + p)/(1 - t)
    # Now we assume the rank == 2.
    m = len(A)
    if style == 'Igusa':
        return (1 - p**-1)*(1 - (m - 1)*p**-1 + (m - 1)*p**-1*t - p**-2*t) / ((1 - p**-1*t)*(1 - p**-2*t**m))
    if style == 'skele':
        return (1 + m*p + (m-1)*p**2 + (m-1 + m*p + p**2)*t)/((1 - t)**2)

# The direct version of the combinatorial skeleton.
def _comb_skele_direct(A, DB=True):
    from sage.all import PolynomialRing, QQ, var
    from Globals import __DEFAULT_t, __DEFAULT_p
    from PosetOps import CharacteristicFunction, PoincarePolynomial, _proper_part
    _, P = CharacteristicFunction(A)
    if DB:
        zeta = _data.get_gen_func(P, 'skele')
        if zeta != None:
            return zeta
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
    zeta = skele/((1 - t)**n)
    if DB:
        _data.save_gen_func(P, 'skele', zeta)
    return zeta

# The direct version of the local Igusa zeta function computation.
def _local_Igusa_direct(A, DB=True):
    from sage.all import PolynomialRing, QQ, var
    from Globals import __DEFAULT_t, __DEFAULT_p
    from PosetOps import CharacteristicFunction, PoincarePolynomial, _proper_part, _Coxeter_poset_data
    p = var(__DEFAULT_p)
    t = var(__DEFAULT_t)
    _, P = CharacteristicFunction(A)
    if DB:
        zeta = _data.get_gen_func(P, 'Igusa')
        if zeta != None:
            return zeta
    # We check to see if we have a type A braid arrangement
    if _Coxeter_poset_data()['A']['hyperplanes'](A.rank()) == len(A):
        if _Coxeter_poset_data()['A']['poset'](A.rank()) == len(P):
            from Constructors import CoxeterArrangement
            B = CoxeterArrangement("A", n=A.rank())
            if P.is_isomorphic(B.intersection_poset()):
                from Braid import BraidArrangementIgusa
                return BraidArrangementIgusa(A.rank())
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
    zeta = skele/reduce(lambda x, y: x*(1 - Y(y)), P._elements[1:], 1)
    if DB:
        _data.save_gen_func(P, 'Igusa', zeta)
    return zeta

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

def _local_Igusa_recurse(A, DB=True):
    from sage.all import PolynomialRing, QQ, var, ZZ
    from Globals import __DEFAULT_t, __DEFAULT_p
    from PosetOps import CharacteristicFunction, PoincarePolynomial, _proper_part, _deletion, _Coxeter_poset_data
    p = var(__DEFAULT_p)
    t = var(__DEFAULT_t)
    char_func, P = CharacteristicFunction(A)
    if DB:
        zeta = _data.get_gen_func(P, 'Igusa')
        if zeta != None:
            return zeta
    # We check to see if we have a type A braid arrangement
    if _Coxeter_poset_data()['A']['hyperplanes'](A.rank()) == len(A):
        if _Coxeter_poset_data()['A']['poset'](A.rank()) == len(P):
            from Constructors import CoxeterArrangement
            B = CoxeterArrangement("A", n=A.rank())
            if P.is_isomorphic(B.intersection_poset()):
                from Braid import BraidArrangementIgusa
                return BraidArrangementIgusa(A.rank())
    P_prop = _proper_part(P)
    hypers = list(filter(lambda x: P.covers('', x), P))
    nHypers = lambda Z: len(list(filter(lambda H: P.le(H, Z), hypers)))
    # p_factor = lambda X: p**len(filter(lambda s: s != '', X.split(' ')))
    t_factor = lambda X: t**nHypers(X)
    elts = P_prop._elements
    factors = map(lambda x: t_factor(x)*char_func(x), elts)
    integrals = map(
        lambda x: LocalIgusaZetaFunction(
            _deletion(A, x, P)[0], method="recursive", database=DB
        ), 
        elts
    )
    zeta = reduce(lambda x, y: x + y[0]*y[1], zip(factors, integrals), 0) + char_func('')
    if A.is_central():
        zeta = p**(-A.dimension())*zeta/(1 - p**(-A.rank())*t**len(A))
    else:
        zeta = p**(-A.dimension())*zeta 
    if DB and A.rank() > 2: 
        _data.save_gen_func(P, 'Igusa', zeta)
    return zeta

def _comb_skele_recurse(A, DB=True):
    from sage.all import PolynomialRing, QQ, var, ZZ
    from Globals import __DEFAULT_t, __DEFAULT_p
    from PosetOps import CharacteristicFunction, PoincarePolynomial, _proper_part, _deletion
    p = var(__DEFAULT_p)
    t = var(__DEFAULT_t)
    char_func, P = CharacteristicFunction(A)
    if DB:
        zeta = _data.get_gen_func(P, 'skele')
        if zeta != None:
            return zeta
    P_prop = _proper_part(P)
    def poincare(x):
        chi = char_func(x)
        d = chi.degree(p)
        return (-p)**d*chi.subs({p: -p**-1})
    elts = P_prop._elements
    factors = map(lambda x: t*poincare(x), elts)
    integrals = map(
        lambda x: CombinatorialSkeleton(
            _deletion(A, x, P)[0], method="recursive", database=DB
        ), 
        elts
    )
    zeta = reduce(lambda x, y: x + y[0]*y[1], zip(factors, integrals), 0) + poincare('')
    if A.is_central():
        zeta = zeta/(1 - t)
    else:
        zeta = zeta 
    if DB and A.rank() > 2: 
        _data.save_gen_func(P, 'skele', zeta)
    return zeta


def CombinatorialSkeleton(A, method="recursive", database=True):
    if A.is_central() and A.rank() <= 2:
        return _small_central(A, 'skele')
    if method == "direct":
        return _comb_skele_direct(A, DB=database)
    else:
        return _comb_skele_recurse(A, DB=database)

def LocalIgusaZetaFunction(A, method="recursive", database=True):
    if A.is_central() and A.rank() <= 2:
        return _small_central(A, 'Igusa')
    if method == "direct":
        return _local_Igusa_direct(A, DB=database)
    else:
        return _local_Igusa_recurse(A, DB=database)

def UniversalGeneratingFunction(A, Map=False):
    return _universal_gen_func(A, Map)
