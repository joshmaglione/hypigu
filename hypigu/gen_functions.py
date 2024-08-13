#
#   Copyright 2021--2024 Joshua Maglione
#
#   Distributed under MIT License
#
from sage.all import QQ, reduce, PolynomialRing, prod, Subsets

from .globals import _PRINT as _print
from .globals import _TIME as _time
from .graded_poset import GradedPoset, _combinatorial_eq_classes

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
    if P.has_top():
        PP = P.subposet(ord_elts[:-1])
    else:
        PP = P
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
    if P.has_top():
        PP = P.subposet(ord_elts[:-1])
    else:
        PP = P
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
    # End recursion
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
    if P.has_top():
        PP = P[1:-1]
    else:
        PP = P[1:]
    fHP = GP.Poincare_polynomial()(Y=Y)
    for x in PP:
        fHP += GP.interval(bottom=x).Poincare_polynomial()(Y=Y) * Ts[x] * _fHP_recursion(GP.interval(top=x), varbs=(Y, Ts))
    if P.has_top():
        fHP = fHP/(1 - Ts[P.top()])
    return fHP

def _cfHP_no_R_label(GP:GradedPoset, numerator=False):
    pass

# Use Theorem 2.7 of Dorpalen-Barry, Maglione, and Stump.
def _cfHP_R_label(GP:GradedPoset, numerator=False):
    Lambda = GP.R_label
    def stat(M, E):
        labels = [Lambda(tup) for tup in zip(M, M[1:])]
        # 'a' == False, 'b' == True
        U = [False] + [tup[0] > tup[1] for tup in zip(labels, labels[1:])]
        stat = 0
        for i, u in enumerate(U):
            if i == 0:
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

def FlagHilbertPoincareSeries(
        arrangement=None,
        matroid=None,
        poset=None,
        verbose=_print
    ):
    GP = GradedPoset(
        arrangement=arrangement,
        matroid=matroid,
        poset=poset
    )
    if verbose:
        print(f"{_time()}Computing the flag Hilbert--Poincaré series")
    return _fHP_recursion(GP)

def CoarseFlagHilbertPoincareSeries(
        arrangement=None,
        matroid=None,
        poset=None,
        R_label=None,
        verbose=_print,
        numerator=False
    ):
    GP = GradedPoset(
        arrangement=arrangement,
        matroid=matroid,
        poset=poset,
        R_label=R_label
    )
    if verbose:
        print(f"{_time()}Computing the coarse flag Hilbert--Poincaré series")
    if R_label is not None or GP.R_label is not None:
        return _cfHP_R_label(GP, numerator=numerator)
    return _cfHP_no_R_label(GP, numerator=numerator)