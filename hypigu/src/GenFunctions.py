#
#   Copyright 2021 Joshua Maglione
#
#   Distributed under MIT License
#

from .Database import internal_database as _data
from .Globals import __PRINT as _print
from .Globals import __TIME as _time
from functools import reduce as _reduce


# A function to return a poincare function.
def _Poincare_polynomial(L, sub=None):
    from sage.all import var
    if sub is None:
        sub = var('Y')

    def poincare(x):
        pi = L.restriction(x).Poincare_polynomial()
        try:
            return pi.subs({pi.variables()[0]: sub})
        except AttributeError:  # In case pi is a constant.
            return pi
    return poincare


# The complete solutions for small central arrangements of rank <= 2.
def _small_central(A, style, numerator=False):
    from sage.all import var

    p = var('q')
    t = var('t')
    Y = var('Y')
    T = var('T')
    if A.rank() == 1:
        if style == 'Igusa':
            return (1 - p**-1)/(1 - p**-1*t)
        if style == 'skele':
            if numerator:
                return 1 + Y
            else:
                return (1 + Y)/(1 - T)
    # Now we assume the rank == 2.
    m = len(A)
    if style == 'Igusa':
        return (1 - p**-1)*(1 - (m - 1)*p**-1 + (m - 1)*p**-1*t - p**-2*t) / ((1 - p**-1*t)*(1 - p**-2*t**m))
    if style == 'skele':
        if numerator:
            return 1 + m*Y + (m-1)*Y**2 + (m-1 + m*Y + Y**2)*T
        else:
            return (1 + m*Y + (m-1)*Y**2 + (m-1 + m*Y + Y**2)*T)/((1 - T)**2)


# The direct version of the universal generating function computation.
def _universal(L, anayltic=False, atom=False):
    from sage.all import var
    from .LatticeFlats import _subposet

    # Set up the potential substitutions for T -- as defined in Maglione--Voll.
    if anayltic:
        q = var('q')
        Y = -q**(-1)
        P = L.poset
        t_name = lambda x: var("t" + str(x))
        if atom:
            atoms = P.upper_covers(P.bottom())

            def T_data(x):
                under_poset = _subposet(P, x, lambda z: P.lower_covers(z))
                elts = filter(lambda y: y in under_poset, atoms)
                ts = map(t_name, elts)
                return _reduce(lambda x, y: x*y, ts, q**(-P.rank(x)))
        else:
            def T_data(x):
                under_poset = _subposet(P, x, lambda z: P.lower_covers(z))
                elts = list(under_poset._elements)
                elts.remove(P.bottom())
                ts = map(t_name, elts)
                return _reduce(lambda x, y: x*y, ts, q**(-P.rank(x)))
    else:
        T_data = lambda x: var("T" + str(x))
        Y = var('Y')

    T = {x: T_data(x) for x in L.poset._elements}

    # Base cases for recursion.
    if L.poset.has_top() and L.poset.rank() == 2:
        elts = L.proper_part_poset()._elements
        merge = lambda x, y: x + (1 + Y)**2*T[y]/(1 - T[y])
        one = L.poset.top()
        return _reduce(merge, elts, (1 + Y)*(1 + (len(elts) - 1)*Y))/(1-T[one])
    if L.poset.rank == 1:
        elts = list(L.poset._elements).remove(L.poset.bottom())
        merge = lambda x, y: x + (1 + Y)*T[y]/(1 - T[y])
        return _reduce(merge, elts, 1 + len(elts)*Y)

    P = L.proper_part_poset()
    poincare = _Poincare_polynomial(L, sub=Y)
    recurse = lambda M: _universal(M, anayltic=anayltic, atom=atom)
    num_dat = lambda x: poincare(x)*T[x]*recurse(L.subarrangement(x))
    factors = map(num_dat, P._elements)
    HP = _reduce(lambda x, y: x + y, factors, poincare(L.poset.bottom()))
    if L.poset.has_top():
        HP = HP/(1 - T[L.poset.top()])
    return HP


def _Igusa_zeta_function(L, DB=True, verbose=_print):
    from sage.all import var
    from .Constructors import CoxeterArrangement
    from .Braid import BraidArrangementIgusa
    from .LatticeFlats import LatticeOfFlats, _Coxeter_poset_data

    P = L.poset
    q = var('q')
    t = var('t')

    # Base cases for recursion.
    if P.has_top() and P.rank() == 2:
        m = len(P) - 2
        return (1 - q**-1)*(1 - (m-1)*q**-1 + m*(1 - q**-1)*q**-1*t/(1 - q**-1*t))/(1 - q**-2*t**m)
    if P.rank() == 1:
        m = len(P) - 1
        return 1 - m*q**-1 + m*(1 - q**-1)*q**-1*t/(1 - q**-1*t)

    if DB:
        zeta = _data.get_gen_func(P, 'Igusa')
        if zeta is not None:
            return zeta
    # We check to see if we have a type A braid arrangement.
    # We can compute these *extremely* quickly.
    if _Coxeter_poset_data()['A']['hyperplanes'](P.rank()) == len(L.atoms()):
        if _Coxeter_poset_data()['A']['poset'](P.rank()) == len(P):
            B = CoxeterArrangement("A" + str(P.rank()))
            if P.is_isomorphic(LatticeOfFlats(B).poset):
                return BraidArrangementIgusa(P.rank())

    poincare = _Poincare_polynomial(L, sub=-q**(-1))
    t_factor = lambda X: t**len(L.flat_labels[X])
    x_factor = lambda x: poincare(x)*t_factor(x)*q**(-P.rank(x))
    eq_elt_data = L._combinatorial_eq_elts()
    factors = map(lambda x: x[1]*x_factor(x[0]), eq_elt_data)
    integrals = map(lambda x: _Igusa_zeta_function(x[2], DB=DB), eq_elt_data)
    pi = poincare(P.bottom())
    zeta = _reduce(lambda x, y: x + y[0]*y[1], zip(factors, integrals), 0) + pi
    if P.has_top():
        zeta = zeta/(1 - q**(-P.rank())*t**len(L.atoms()))
    if DB and P.rank() > 2:
        _data.save_gen_func(P, 'Igusa', zeta)
    return zeta


def _top_zeta_function_uni(L, DB=True, verbose=_print):
    from sage.all import var

    P = L.poset
    s = var('s')
    C = 1*L.poset.has_top()

    # Base cases for recursion.
    if P.has_top() and P.rank() == 2:
        m = len(P) - 2
        return (2 + (2 - m)*s)/((2 + m*s)*(1 + s))
    if P.rank() == 1:
        m = len(P) - 1
        return (1 + (1 - m)*s)/(1 + s)

    poincare = _Poincare_polynomial(L)
    Y = poincare(P.bottom()).variables()[0]
    pi_circ = lambda x: (poincare(x)/(1 + Y)**C).factor().simplify().subs({Y: -1})
    eq_elt_data = L._combinatorial_eq_elts()
    factors = map(lambda x: x[1]*pi_circ(x[0]), eq_elt_data)
    integrals = map(lambda x: _top_zeta_function_uni(x[2], DB=DB), eq_elt_data)
    pi = pi_circ(P.bottom())
    zeta = _reduce(lambda x, y: x + y[0]*y[1], zip(factors, integrals), 0) + pi
    if C == 1:
        zeta = zeta/(P.rank() + len(L.atoms())*s)

    return zeta


def _top_zeta_function_mul(L, DB=True, verbose=_print, atom=False):
    from sage.all import var
    from .LatticeFlats import _subposet

    P = L.poset
    C = 1 * L.poset.has_top()

    s_name = lambda x: var("s" + str(x))
    if atom:
        if atom:
            atoms = P.upper_covers(P.bottom())

            def s_data(x):
                under_poset = _subposet(P, x, lambda z: P.lower_covers(z))
                elts = filter(lambda y: y in under_poset, atoms)
                ts = map(s_name, elts)
                return _reduce(lambda x, y: x + y, ts, 0)
    else:
        def s_data(x):
            under_poset = _subposet(P, x, lambda z: P.lower_covers(z))
            elts = list(under_poset._elements)
            elts.remove(P.bottom())
            ts = map(s_name, elts)
            return _reduce(lambda x, y: x + y, ts, 0)

    S = {x: s_data(x) for x in P._elements}

    # Base cases for recursion.
    add_em = lambda x, y: x + y
    if P.has_top() and P.rank() == 2:
        atms = P.upper_covers(P.bottom())
        m = len(atms)
        elt_dat = lambda x: 1/(1 + S[x])
        return _reduce(add_em, map(elt_dat, atms), 2 - m)/(2 + S[P.top()])
    if P.rank() == 1:
        atms = P.upper_covers(P.bottom())
        m = len(atms)
        elt_dat = lambda x: 1/(1 + S[x])
        return _reduce(add_em, map(elt_dat, atms), 1 - m)

    poincare = _Poincare_polynomial(L)
    Y = poincare(P.bottom()).variables()[0]
    pi_circ = lambda x: (poincare(x)/(1 + Y)**C).factor().simplify().subs({Y: -1})
    x_factor = lambda x: pi_circ(x)
    prop_elts = L.proper_part_poset()._elements
    factors = map(lambda x: x_factor(x), prop_elts)
    integrals = map(lambda x: _top_zeta_function_mul(L.subarrangement(x), DB=DB, atom=atom), prop_elts)
    pi = pi_circ(P.bottom())
    zeta = _reduce(lambda x, y: x + y[0]*y[1], zip(factors, integrals), 0) + pi
    if P.has_top():
        zeta = zeta/(P.rank() + S[P.top()])

    return zeta


def _comb_skele(L, DB=True, verbose=_print):
    from sage.all import var
    P = L.poset
    Y = var('Y')
    T = var('T')

    if P.has_top():
        if P.rank() == 1:
            return (1 + Y)/(1 - T)
        if P.rank() == 2:
            m = len(P) - 2
            return (1 + m*Y + (m - 1)*Y**2 + (m - 1 + m*Y + Y**2)*T)/(1 - T)**2
    if DB:
        if verbose:
            print(_time() + "Checking database.")
        zeta = _data.get_gen_func(P, 'skele')
        if zeta is not None:
            return zeta
        if verbose:
            print("\tDone.")

    poincare = _Poincare_polynomial(L)
    if verbose:
        print(_time() + "Gleaning structure from poset.")
    eq_elt_data = L._combinatorial_eq_elts()
    if verbose:
        print("\tDone.")
        print(_time() + "Lattice points: {0},  Relevant points: {1}".format(len(P), len(eq_elt_data)))
    factors = map(lambda x: x[1]*T*poincare(x[0]), eq_elt_data)
    if verbose:
        print(_time() + "Recursing...")
    integrals = map(lambda x: _comb_skele(x[2], DB=DB), eq_elt_data)
    if verbose:
        print(_time() + "Putting everything together...")
    pi = poincare(P.bottom())
    zeta = _reduce(lambda x, y: x + y[0]*y[1], zip(factors, integrals), 0) + pi
    if P.has_top():
        zeta = zeta/(1 - T)
    if DB and P.rank() > 2:
        _data.save_gen_func(P, 'skele', zeta)

    return zeta


# Given a polynomial, return a hyperplane arrangement equivalent to the linear
# factors of f.
def _parse_poly(f):
    from sage.all import SR, QQ, HyperplaneArrangements, Matrix
    if type(f) == str:
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
    is_lin = lambda g: all(map(lambda x: g.degree(x) <= 1, g.variables()))
    if not all(map(is_lin, F)):
        raise ValueError("Expected product of linear factors.")

    varbs = f.variables()
    varbs_str = tuple(map(lambda x: str(x), varbs))
    HH = HyperplaneArrangements(K, varbs_str)

    def poly_vec(g):
        c = K(g.subs({x: 0 for x in g.variables()}))
        return tuple([c] + [K(g.coefficient(x)) for x in varbs])

    F_vec = tuple(map(poly_vec, F))
    A = HH(Matrix(K, F_vec))

    # This scrambles the hyperplanes, so we need to scramble M in the same way.
    A_vec = tuple(map(lambda H: tuple(H.coefficients()), A.hyperplanes()))
    perm = tuple([F_vec.index(v) for v in A_vec])
    M_new = tuple([M[i] for i in perm])

    return A, M_new


def CoarseFlagHPSeries(A=None, lattice_of_flats=None, int_poset=None, matroid=None, numerator=False, verbose=_print):
    from .LatticeFlats import LatticeOfFlats

    if matroid is None:
        try:
            if A.is_central() and A.rank() <= 2:
                return _small_central(A, 'skele', numerator=numerator)
        except AttributeError:
            raise TypeError("object is not a hyperplane arrangement.")
    if lattice_of_flats is None:
        if verbose:
            print("{0}Building lattice of flats".format(_time()))
        if matroid is None:
            L = LatticeOfFlats(A, poset=int_poset)
        else:
            L = LatticeOfFlats(matroid=matroid)
    else:
        L = lattice_of_flats

    if verbose:
        print("{0}Computing coarse flag Hilbert--Poincare series".format(_time()))
    cfHP = _comb_skele(L)

    if numerator:
        D = cfHP.numerator_denominator()[1]
        T = D.variables()[0]
        if D == (T - 1)**(L.poset.rank()):
            e = -1
        if D == (1 - T)**(L.poset.rank()):
            e = 1
        return e*(cfHP*D).factor()
    else:
        return cfHP


def IgusaZetaFunction(X=None, lattice_of_flats=None, int_poset=None, matroid=None, verbose=_print):
    from .LatticeFlats import LatticeOfFlats
    from sage.all import var

    HPA = True
    if matroid is None:
        try:
            # Check if a hyperplane arrangement.
            _ = X.hyperplanes()
            A = X
        except AttributeError:
            # Not an HPA; deal with polynomial input.
            A, M = _parse_poly(X)
            if verbose:
                print("{0}Constructed a hyperplane arrangement".format(_time()))
            HPA = False

    if lattice_of_flats is None:
        if verbose:
            print("{0}Building lattice of flats".format(_time()))
        if matroid is None:
            L = LatticeOfFlats(A, poset=int_poset)
        else:
            L = LatticeOfFlats(matroid=matroid)
    else:
        L = lattice_of_flats

    if not HPA:
        if list(M) == [1]*len(M):
            if verbose:
                print("{0}Computing Igusa's zeta function".format(_time()))
            return _Igusa_zeta_function(L)
        else:
            if verbose:
                print("{0}Computing the atom zeta function".format(_time()))
            Z = _universal(L, anayltic=True, atom=True)
            t = var('t')
            SUB = {var('t' + str(k+1)): t**(M[k]) for k in range(len(M))}
            return Z.subs(SUB)

    if verbose:
        print("{0}Computing Igusa's zeta function".format(_time()))
    return _Igusa_zeta_function(L)


def TopologicalZetaFunction(X=None, lattice_of_flats=None, int_poset=None, verbose=_print, multivariate=False, atom=False, matroid=None):
    from .LatticeFlats import LatticeOfFlats
    from sage.all import var

    HPA = True
    if matroid is None:
        try:
            # Check if a hyperplane arrangement.
            _ = X.hyperplanes()
            A = X
        except AttributeError:
            # Not an HPA; deal with polynomial input.
            A, M = _parse_poly(X)
            if verbose:
                print("{0}Constructed a hyperplane arrangement".format(_time()))
            HPA = False

    if lattice_of_flats is None:
        if matroid is None:
            if verbose:
                print("{0}Building lattice of flats".format(_time()))
            L = LatticeOfFlats(A, poset=int_poset)
        else:
            L = LatticeOfFlats(matroid=matroid)
    else:
        L = lattice_of_flats

    if verbose:
        print("{0}Computing the topological zeta function".format(_time()))

    if not HPA:
        if list(M) == [1] * len(M):
            return _top_zeta_function_uni(L)
        else:
            Z = _top_zeta_function_mul(L, atom=True)
            s = var('s')
            SUB = {var('s' + str(k+1)): M[k]*s for k in range(len(M))}
            return Z.subs(SUB)

    if not multivariate:
        return _top_zeta_function_uni(L)

    return _top_zeta_function_mul(L, atom=atom)


def AnalyticZetaFunction(A=None, lattice_of_flats=None, int_poset=None, matroid=None, verbose=_print):
    from .LatticeFlats import LatticeOfFlats

    if lattice_of_flats is None:
        if matroid is None:
            if verbose:
                print("{0}Building lattice of flats".format(_time()))
            L = LatticeOfFlats(A, poset=int_poset)
        else:
            L = LatticeOfFlats(matroid=matroid)
    else:
        L = lattice_of_flats

    if verbose:
        print("{0}Computing the analytic zeta function".format(_time()))
    return _universal(L, anayltic=True)


def AtomZetaFunction(A=None, lattice_of_flats=None, int_poset=None, matroid=None, verbose=_print):
    from .LatticeFlats import LatticeOfFlats

    if lattice_of_flats is None:
        if matroid is None:
            if verbose:
                print("{0}Building lattice of flats".format(_time()))
            L = LatticeOfFlats(A, poset=int_poset)
        else:
            L = LatticeOfFlats(matroid=matroid)
    else:
        L = lattice_of_flats

    if verbose:
        print("{0}Computing the atom zeta function".format(_time()))
    return _universal(L, anayltic=True, atom=True)


def FlagHilbertPoincareSeries(A=None, lattice_of_flats=None, int_poset=None, matroid=None, verbose=_print):
    from .LatticeFlats import LatticeOfFlats

    if lattice_of_flats is None:
        if matroid is None:
            if verbose:
                print("{0}Building lattice of flats".format(_time()))
            L = LatticeOfFlats(A, poset=int_poset)
        else:
            L = LatticeOfFlats(matroid=matroid)
    else:
        L = lattice_of_flats

    if verbose:
        print("{0}Computing the flag Hilbert--Poincare series".format(_time()))
    return _universal(L)
