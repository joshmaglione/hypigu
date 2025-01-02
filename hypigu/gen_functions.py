#
#   Copyright 2021--2024 Joshua Maglione
#
#   Distributed under MIT License
#
from sage.rings.rational_field import Q as QQ
from _functools import reduce
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.misc.misc_c import prod
from sage.combinat.subset import Subsets
from sage.rings.integer_ring import Z as ZZ
from sage.rings.polynomial.polynomial_ring import polygens

from .globals import my_print, verbose
from .graded_poset import GradedPoset, proper_part

def combinatorial_eq_classes(GP:GradedPoset):
    def classes(P, elts):
        top = lambda x: P.subposet(P.principal_order_filter(x))
        bot = lambda x: P.subposet(P.principal_order_ideal(x))
        cl = []
        for x in elts:
            matched = False
            for c in cl:
                if P.rank_function()(x) != P.rank_function()(c[0]):
                    continue
                if len(P.upper_covers(x)) != len(P.upper_covers(c[0])):
                    continue
                if len(P.lower_covers(x)) != len(P.lower_covers(c[0])):
                    continue
                if top(x).is_isomorphic(top(c[0])) and bot(x).is_isomorphic(bot(c[0])):
                    c.append(x)
                    matched = True
                    break
            if not matched:
                cl.append([x])
        return cl
    return classes(GP.poset, proper_part(GP.poset))

# Construct fHP directly from the definition. 
def _fHP_chains(GP:GradedPoset):
    P = GP.poset
    ord_elts = [list(filter(
        lambda x: P.rank_function()(x) == r, P
    )) for r in range(P.rank() + 1)]
    ord_elts = reduce(lambda x, y: x + y, ord_elts, [])
    ord_elts = ord_elts[1:]
    Ts = [f"T{i + 1}" for i, _ in enumerate(ord_elts)]
    PR = PolynomialRing(QQ, ["Y"] + Ts)
    Y = PR.gens()[0]
    Ts = {x: PR.gens()[i + 1] for i, x in enumerate(ord_elts)}
    geoprog = lambda S: Ts[S]/(1 - Ts[S])
    def chain_Poincare(F):
        print(F)
        if len(F) == 0:
            return GP.Poincare_polynomial()(Y=Y)
        tups = [(None, F[0])] + list(zip(F, F[1:])) + [(F[-1], None)]
        Poincares = [GP.interval(*t).Poincare_polynomial()(Y=Y) for t in tups]
        return reduce(lambda x, y: x*y, Poincares, PR.one())
    PP = proper_part(P)
    terms = [chain_Poincare(F) * reduce(
        lambda x, y: x*y, map(geoprog, F), PR.one()
    ) for F in PP.chains()]
    fHP = sum(terms)
    if P.has_top():
        fHP = fHP/(1 - Ts[P.top()])
    return fHP


# Construct the numerator of fHP directly from the definition.
def _fHP_numerator(GP:GradedPoset):
    P = GP.poset
    ord_elts = [list(filter(
        lambda x: P.rank_function()(x) == r, P
    )) for r in range(P.rank() + 1)]
    ord_elts = reduce(lambda x, y: x + y, ord_elts, [])
    ord_elts = ord_elts[1:]
    Ts = [f"T{i + 1}" for i, _ in enumerate(ord_elts)]
    PR = PolynomialRing(QQ, ["Y"] + Ts)
    Y = PR.gens()[0]
    Ts = {x: PR.gens()[i + 1] for i, x in enumerate(ord_elts)}
    def chain_Poincare(F):
        print(F)
        if len(F) == 0:
            return GP.Poincare_polynomial()(Y=Y)
        tups = [(None, F[0])] + list(zip(F, F[1:])) + [(F[-1], None)]
        Poincares = [GP.interval(*t).Poincare_polynomial()(Y=Y) for t in tups]
        return reduce(lambda x, y: x*y, Poincares, PR.one())
    PP = proper_part(P)
    PP_set = frozenset(PP)
    numer = lambda F: (-1)**(len(F) % 2) * prod(
        Ts[S] for S in F
    ) * prod(
        Ts[S] - 1 for S in PP_set.difference(frozenset(F))
    )
    numerator = sum(chain_Poincare(F) * numer(F) for F in PP.chains())
    fHP = numerator / prod(1 - T for T in Ts.values())
    return fHP

# Based off an argument very similar to Lemma 2.2 of Maglione--Voll.
def _fHP_recursion(GP:GradedPoset, varbs=None):
    # Initial Setup
    P = GP.poset
    if varbs is None:
        ord_elts = [list(filter(
            lambda x: P.rank_function()(x) == r, P
        )) for r in range(P.rank() + 1)]
        ord_elts = reduce(lambda x, y: x + y, ord_elts, [])
        ord_elts = ord_elts[1:]
        Ts = [f"T{i + 1}" for i, _ in enumerate(ord_elts)]
        PR = PolynomialRing(QQ, ["Y"] + Ts)
        Y = PR.gens()[0]
        Ts = {x: PR.gens()[i + 1] for i, x in enumerate(ord_elts)}
    else: 
        Y, Ts = varbs
    # Base cases
    if P.rank() == 1: # Proposition 4.1 of Maglione--Voll.
        m = len(GP.atoms())
        return 1 + m*Y + (1 + Y)*sum(Ts[S]/(1 - Ts[S]) for S in GP.atoms())
    if P.rank() == 2 and P.has_top(): # Proposition 4.2 of Maglione--Voll.
        m = len(GP.atoms())
        return (1 + Y)/(1 - Ts[P.top()]) * (1 + (m - 1)*Y + (1 + Y)*sum(
            Ts[S]/(1 - Ts[S]) for S in GP.atoms())
        )
    # Determine if seen before
    # 
    # Recursion
    PP = proper_part(P)
    fHP = GP.Poincare_polynomial()(Y=Y)
    for x in PP:
        fHP += GP.interval(bottom=x).Poincare_polynomial()(Y=Y) * Ts[x] * _fHP_recursion(GP.interval(top=x), varbs=(Y, Ts))
    if P.has_top():
        fHP = fHP/(1 - Ts[P.top()])
    return fHP

def _cfHP_no_R_label(GP:GradedPoset, varbs=None):
    # Initial Setup
    P = GP.poset
    if varbs is None:
        PR = PolynomialRing(QQ, 'Y,T')
        Y, T = PR.gens()
    else:
        Y, T = varbs
    # Base cases
    if P.rank() == 1: # Proposition 4.1 of Maglione--Voll.
        m = len(GP.atoms())
        return (1 + m*Y + (m - 1)*T) / (1 - T)
    if P.rank() == 2 and P.has_top(): # Proposition 4.2 of Maglione--Voll.
        m = len(GP.atoms())
        return (1 + m*Y + (m - 1)*Y**2 + (m - 1 + m*Y + Y**2)*T) / (1 - T)**2
    # Determine if seen before
    # 
    # Recursion
    Classes = combinatorial_eq_classes(GP)
    fHP = GP.Poincare_polynomial()(Y=Y)
    for C in Classes:
        x = C[0]
        M = len(C)
        fHP += M * GP.interval(bottom=x).Poincare_polynomial()(Y=Y) * T * _cfHP_no_R_label(GP.interval(top=x), varbs=(Y, T))
    if P.has_top():
        fHP = fHP/(1 - T)
    return fHP

# Use Theorem 2.7 & Corollary 2.22 of Dorpalen-Barry, Maglione, and Stump.
def _cfHP_R_label(GP:GradedPoset, numerator=False):
    Lambda = GP.R_label
    def stat(M, E):
        labels = [Lambda(tup) for tup in zip(M, M[1:])]
        # 'a' == False, 'b' == True
        U = [False] + [tup[0] > tup[1] for tup in zip(labels, labels[1:])]
        stat = 0
        for i, u in enumerate(U):
            if i == 0: # iota substitution
                continue
            if (not u and i in E) or (u and i - 1 not in E):
                stat += 1
        return stat
    P = GP.poset
    PR = PolynomialRing(QQ, 'Y,T')
    Y, T = PR.gens()
    r = P.rank()
    numer = sum(
        Y**(len(E)) * T**(stat(M, E)) for E in Subsets(range(r)) for M in P.maximal_chains()
    )
    if numerator:
        return numer
    return numer / (1 - T)**r

def _IZF_recursion(GP:GradedPoset, varbs=None):
    # Initial Setup
    P = GP.poset
    # y = q^{-1} and t = q^{-s}
    if varbs is None:
        PR = PolynomialRing(QQ, 'y, t')
        y, t = PR.gens()
    else:
        y, t = varbs
    # Base cases
    if P.rank() == 1: # Proposition 4.1 of Maglione--Voll.
        m = len(GP.atoms())
        return (1 - m*y + (m - 1)*y*t) / (1 - y*t)
    if P.rank() == 2 and P.has_top(): # Proposition 4.2 of Maglione--Voll.
        m = len(GP.atoms())
        return (1 - m*y + (m - 1)*y**2 + (m - 1 - m*y + y**2)*y*t) / ((1 - y*t) * (1 - y**2*t**m))
    # Determine if seen before
    # 
    # Recursion
    Classes = combinatorial_eq_classes(GP)
    Z = GP.Poincare_polynomial()(Y=-y)
    for C in Classes:
        x = C[0]
        M = len(C)
        l = len([a for a in GP.atoms() if P.le(a, x)])
        r = P.rank_function()(x)
        Z += M * GP.interval(bottom=x).Poincare_polynomial()(Y=-y) * y**r * t**l * _IZF_recursion(GP.interval(top=x), varbs=(y, t))
    if P.has_top():
        R = P.rank()
        A = len(GP.atoms())
        Z = Z/(1 - y**R * t**A)
    return Z

def _TZF_chains(GP:GradedPoset):
    P = GP.poset
    s = polygens(ZZ, 's')[0]
    def chain_Poincare_circ(F):
        if len(F) == 0:
            pi = GP.Poincare_polynomial()
        else:
            tups = [(None, F[0])] + list(zip(F, F[1:])) + [(F[-1], None)]
            Poincares = [GP.interval(*t).Poincare_polynomial() for t in tups]
            pi = reduce(lambda x, y: x*y, Poincares, 1)
        Y = pi.parent().gens()[0]
        pi = pi/(1 + Y)**(len(F) + int(P.has_top()))
        return pi(Y=-1)
    gx = lambda x: P.rank_function()(x) + len([a for a in GP.atoms() if P.le(a, x)])*s
    PP = proper_part(P)
    Z = sum(
        chain_Poincare_circ(F) * prod(1/gx(S) for S in F) for F in PP.chains()
    )
    if P.has_top():
        Z = Z/(P.rank() + len(GP.atoms())*s)
    return Z

def _TZF_recursion(GP:GradedPoset, varbs=None):
    # Initial Setup
    P = GP.poset
    if varbs is None:
        PR = PolynomialRing(ZZ, 's')
        s = PR.gens()[0]
    else:
        s = varbs[0]
    def Poin_circ(x):
        pi = GP.interval(bottom=x).Poincare_polynomial()
        if P.has_top():
            Y = pi.parent().gens()[0]
            return pi/(1 + Y)
        return pi
    # Base cases
    if P.rank() == 1: # Proposition 4.1 of Maglione--Voll.
        m = len(GP.atoms())
        return (1 - (m - 1)*s) / (1 + s)
    if P.rank() == 2 and P.has_top(): # Proposition 4.2 of Maglione--Voll.
        m = len(GP.atoms())
        return (2 - (m - 2)*s) / ((1 + s) * (2 + m*s))
    # Determine if seen before
    # 
    # Recursion
    Classes = combinatorial_eq_classes(GP)
    Z = Poin_circ(P.bottom())(Y=-1)
    for C in Classes:
        x = C[0]
        M = len(C)
        Z += M * Poin_circ(x)(Y=-1) * _TZF_recursion(GP.interval(top=x), varbs=(s,))
    if P.has_top():
        R = P.rank()
        A = len(GP.atoms())
        Z = Z/(R + A*s)
    return Z


def FlagHilbertPoincareSeries(
    arrangement=None,
    matroid=None,
    poset=None,
    verbose=verbose
):
    GP = GradedPoset(arrangement=arrangement, matroid=matroid, poset=poset)
    my_print(verbose, "Computing the flag Hilbert--Poincaré series.")
    return _fHP_recursion(GP)

def CoarseFHPSeries(
    arrangement=None,
    matroid=None,
    poset=None,
    R_label=None,
    verbose=verbose,
    numerator=False,
    method='recursion',
):
    GP = GradedPoset(
        arrangement=arrangement,
        matroid=matroid,
        poset=poset,
        R_label=R_label
    )
    my_print(verbose, "Computing the coarse flag Hilbert--Poincaré series.")
    if (R_label is not None or GP.R_label is not None) and method.lower() == 'r-label':
        my_print(verbose, "Using the R-label.", level=1)
        return _cfHP_R_label(GP, numerator=numerator)
    cfHP = _cfHP_no_R_label(GP)
    if numerator:
        _, T = cfHP.parent().gens()
        cfHP = cfHP.numerator()*((1 - T)**GP.poset.rank() / cfHP.denominator())
    return cfHP

def IgusaZetaFunction(
    arrangement=None,
    matroid=None,
    poset=None,
    verbose=verbose
):
    GP = GradedPoset(arrangement=arrangement, matroid=matroid, poset=poset)
    my_print(verbose, "Computing the Igusa zeta function.")
    return _IZF_recursion(GP)

def TopologicalZetaFunction(
    arrangement=None,
    matroid=None,
    poset=None,
    verbose=verbose
):
    GP = GradedPoset(arrangement=arrangement, matroid=matroid, poset=poset)
    my_print(verbose, "Computing the Igusa zeta function.")
    return _TZF_recursion(GP)