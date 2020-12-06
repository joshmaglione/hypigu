#
#   Copyright 2020 Joshua Maglione 
#
#   Distributed under MIT License
#

from .Database import internal_database as _data
from .Globals import __PRINT as _print
from .Globals import __TIME as _time
from .LatticeFlats import LatticeOfFlats as _LOF
from functools import reduce as _reduce


# The complete solutions for small central arrangements of rank <= 2.
def _small_central(A, style):
    from .Globals import __DEFAULT_p, __DEFAULT_t 
    from sage.all import var

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
    from .Globals import __DEFAULT_t, __DEFAULT_p
    from .PosetOps import CharacteristicFunction, PoincarePolynomial, _proper_part
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
        Z_in = _reduce(lambda x, y: x*Y(y), C, 1)
        C_comp = filter(lambda z: not z in C, P_prop._elements)
        Z_out = _reduce(lambda x, y: x*(1 - Y(y)), C_comp, 1)
        skele += pi.subs({p: Yvar})*Z_in*Z_out
    Skeleton = skele/_reduce(lambda x, y: x*(1 - Y(y)), P._elements[1:], 1)
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

def _local_Igusa(A, DB=True, poset=None, OG=None):
    from sage.all import PolynomialRing, QQ, var, ZZ
    from .Globals import __DEFAULT_t, __DEFAULT_p
    from .PosetOps import CharacteristicFunction, PoincarePolynomial, _deletion, _Coxeter_poset_data, IntersectionPoset, _equiv_elts
    from .Constructors import CoxeterArrangement
    from .Braid import BraidArrangementIgusa

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
        lambda x: _local_Igusa(
            _deletion(A, x[0], P, poset=False, OG=OG), 
            DB=DB,
            poset=x[2],
            OG=OG
        ), 
        eq_elt_data
    )
    zeta = _reduce(lambda x, y: x + y[0]*y[1], zip(factors, integrals), 0) + char_func('')
    if A.is_central():
        zeta = p**(-A.dimension())*zeta/(1 - p**(-A.rank())*t**len(A))
    else:
        zeta = p**(-A.dimension())*zeta 
    if DB and A.rank() > 2: 
        _data.save_gen_func(P, 'Igusa', zeta)
    return zeta


def _comb_skele(L, DB=True, verbose=_print):
    from sage.all import PolynomialRing, QQ, var, ZZ
    from .Globals import __DEFAULT_t, __DEFAULT_p

    P = L.poset
    p = var(__DEFAULT_p)
    t = var(__DEFAULT_t)

    if P.has_top():
        if P.rank() == 1:
            return (1 + p)/(1 - t)
        if P.rank() == 2:
            m = len(P) - 2
            return (1 + m*p + (m - 1)*p**2 + (m - 1 + m*p + p**2)*t)/(1 - t)**2
    if DB:
        if verbose:
            print(_time() + "Checking database.")
        zeta = _data.get_gen_func(P, 'skele')
        if zeta != None:
            return zeta
        if verbose:
            print("\tDone.")

    def poincare(x):
        pi = L.restriction(x).Poincare_polynomial()
        return pi.subs({pi.variables()[0] : p})
    if verbose: 
        print(_time() + "Gleaning structure from poset.")
    eq_elt_data = L._combinatorial_eq_elts()
    if verbose:
        print("\tDone.")
        print(_time() + "Found the following basic structure:\n%s" % (eq_elt_data))
    factors = map(lambda x: x[1]*t*poincare(x[0]), eq_elt_data)
    if verbose:
        print(_time() + "Recursing...")
    integrals = map(lambda x: _comb_skele(x[2], DB=DB), eq_elt_data)
    if verbose:
        print(_time() + "Putting everything together...")
    pi = L.Poincare_polynomial()
    pi = pi.subs({pi.variables()[0] : p})
    zeta = _reduce(lambda x, y: x + y[0]*y[1], zip(factors, integrals), 0) + pi
    if P.has_top():
        zeta = zeta/(1 - t)
    else:
        zeta = zeta 
    if DB and P.rank() > 2: 
        _data.save_gen_func(P, 'skele', zeta)
        
    return zeta




def CombinatorialSkeleton(A, database=True, lattice_of_flats=None, int_poset=None, verbose=_print):
    if A.is_central() and A.rank() <= 2:
        return _small_central(A, 'skele')
    if lattice_of_flats == None:
        L = _LOF(A, poset=int_poset)
    else:
        L = lattice_of_flats

    return _comb_skele(L, DB=database)

def LocalIgusaZetaFunction(A, database=True, int_poset=None):
    if A.is_central() and A.rank() <= 2:
        return _small_central(A, 'Igusa')
    return _local_Igusa(A, DB=database, poset=int_poset, OG=A)

def UniversalGeneratingFunction(A, Map=False, int_poset=None):
    return _universal_gen_func(A, Map)
