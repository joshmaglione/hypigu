#
#   Copyright 2020 Joshua Maglione 
#
#   Distributed under MIT License
#

from Database import internal_database as _data
from Globals import __PRINT as _print

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

def _local_Igusa_BEST(A, DB=True, poset=None, OG=None):
    from sage.all import PolynomialRing, QQ, var, ZZ
    from Globals import __DEFAULT_t, __DEFAULT_p
    from PosetOps import CharacteristicFunction, PoincarePolynomial, _deletion, _Coxeter_poset_data, IntersectionPoset, _equiv_elts
    from Constructors import CoxeterArrangement
    from Braid import BraidArrangementIgusa

    if A.is_central() and A.rank() <= 2:
        return _small_central(A, 'Igusa')
    p = var(__DEFAULT_p)
    t = var(__DEFAULT_t)
    if poset:
        P = poset
    else:
        P = IntersectionPoset(A)
    if DB:
        zeta = _data.get_gen_func(P, 'Igusa')
        if zeta != None:
            return zeta
    # We check to see if we have a type A braid arrangement. 
    if _Coxeter_poset_data()['A']['hyperplanes'](A.rank()) == len(A):
        if _Coxeter_poset_data()['A']['poset'](A.rank()) == len(P):
            B = CoxeterArrangement("A", n=A.rank())
            if P.is_isomorphic(B.intersection_poset()):
                return BraidArrangementIgusa(A.rank())
    char_func = CharacteristicFunction(A, poset=P)
    def t_factor(X):
        if X == '':
            return 0
        else:
            return len(X.split(' '))
    eq_elt_data = _equiv_elts(P)
    factors = map(lambda x: x[1]*t_factor(x[0])*char_func(x[0]), eq_elt_data)
    integrals = map(
        lambda x: _local_Igusa_BEST(
            _deletion(A, x[0], P, poset=False, OG=OG), 
            DB=DB,
            poset=x[2],
            OG=OG
        ), 
        eq_elt_data
    )
    zeta = reduce(lambda x, y: x + y[0]*y[1], zip(factors, integrals), 0) + char_func('')
    if A.is_central():
        zeta = p**(-A.dimension())*zeta/(1 - p**(-A.rank())*t**len(A))
    else:
        zeta = p**(-A.dimension())*zeta 
    if DB and A.rank() > 2: 
        _data.save_gen_func(P, 'Igusa', zeta)
    return zeta

def _comb_skele_BEST(A, DB=True, poset=None, OG=None, verbose=_print):
    from sage.all import PolynomialRing, QQ, var, ZZ
    from Globals import __DEFAULT_t, __DEFAULT_p
    from PosetOps import CharacteristicFunction, _deletion, IntersectionPoset, _equiv_elts

    if A.is_central() and A.rank() <= 2:
        return _small_central(A, 'skele')
    p = var(__DEFAULT_p)
    t = var(__DEFAULT_t)
    if poset:
        P = poset
    else:
        if verbose: 
            print("Constructing intersection poset.")
        P = IntersectionPoset(A)
        if verbose:
            print("\tDone.")
    if DB:
        if verbose:
            print("Checking database.")
        zeta = _data.get_gen_func(P, 'skele')
        if zeta != None:
            return zeta
        if verbose:
            print("\tDone.")
    char_func = CharacteristicFunction(A, poset=P)
    def poincare(x):
        chi = char_func(x)
        d = chi.degree(p)
        return (-p)**d*chi.subs({p: -p**-1})
    if verbose: 
        print("Gleaning structure from poset.")
    eq_elt_data = _equiv_elts(P)
    if verbose:
        print("\tDone.")
        print("Found the following basic structure:\n%s" % (eq_elt_data))
    factors = map(lambda x: x[1]*t*poincare(x[0]), eq_elt_data)
    if verbose:
        print("Recursing...")
    integrals = map(
        lambda x: _comb_skele_BEST(
            _deletion(A, x[0], P, poset=False, OG=OG), 
            DB=DB,
            poset=x[2],
            OG=OG
        ), 
        eq_elt_data
    )
    if verbose:
        print("Putting everything together...")
    zeta = reduce(lambda x, y: x + y[0]*y[1], zip(factors, integrals), 0) + poincare('')
    if A.is_central():
        zeta = zeta/(1 - t)
    else:
        zeta = zeta 
    if DB and A.rank() > 2: 
        _data.save_gen_func(P, 'skele', zeta)
    return zeta



def CombinatorialSkeleton(A, database=True, int_poset=None):
    return _comb_skele_BEST(A, DB=database, poset=int_poset, OG=A)

def LocalIgusaZetaFunction(A, database=True, int_poset=None):
    return _local_Igusa_BEST(A, DB=database, poset=int_poset, OG=A)

def UniversalGeneratingFunction(A, Map=False, int_poset=None):
    return _universal_gen_func(A, Map)
