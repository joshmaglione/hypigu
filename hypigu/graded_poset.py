#
#   Copyright 2020--2024 Joshua Maglione
#
#   Distributed under MIT License
#
from sage.all import Matrix, Poset, Matroid, QQ, PolynomialRing, DiGraph
from sage.misc.cachefunc import cached_method


class GradedPoset():

    def __init__(
            self,
            arrangement=None,
            matroid=None,
            poset=None,
            R_label=None
        ):
        self.arrangement = arrangement
        self.matroid = matroid
        self.poset = poset
        self.R_label = R_label
        if arrangement is None and matroid is None and poset is None:
            raise ValueError("Expected either an arrangement, matroid, or poset.")
        if poset is not None:
            if not poset.has_bottom():
                raise ValueError("Expected a unique minimal element in poset.")
            if not poset.is_graded(): 
                raise ValueError("Expected a graded poset.")
        if matroid is not None:
            self.poset = matroid.lattice_of_flats()
            self.R_label = lambda e: min(e[1].difference(e[0]))
        if arrangement is not None:
            M = Matrix([H.normal() for H in arrangement])
            self.R_label = lambda e: min(e[1].difference(e[0]))
            if arrangement.is_central(): 
                self.matroid = Matroid(M.transpose())
                self.poset = self.matroid.lattice_of_flats()
            else:
                A_cone = arrangement.cone()
                M_cone = Matrix([H.normal() for H in A_cone])
                i = ([list(x) for x in M_cone]
                     .index([1] + [0]*(M_cone.ncols() - 1)))
                MT = Matroid(M_cone.transpose())
                L = MT.lattice_of_flats()
                P = L.subposet(list(filter(lambda x: not i in x, L)))
                bij = [None for _ in range(len(A_cone))]
                for j, H in enumerate(A_cone):
                    if i != j:
                        bij[j] = list(M).index(H.normal()[1:])  
                    else:
                        bij[j] = None
                set_bij = lambda S: frozenset([bij[s] for s in S])
                E = P.cover_relations_graph().edges()
                self.poset = Poset(DiGraph(
                    [(set_bij(e[0]), set_bij(e[1])) for e in E]
                ))

    def __repr__(self):
        if self.R_label:
            start = f"An R-labeled graded poset with {len(self.poset)} elements"
        else:
            start = f"A graded poset with {len(self.poset)} elements"
        if self.arrangement:
            return f"{start} built from the intersection poset of\n{self.arrangement}"
        if self.matroid:
            return f"{start} built from the lattice of flats of\n{self.matroid}"
        return start

    def atoms(self):
        return self.poset.upper_covers(self.poset.bottom())

    def interval(self, bottom=None, top=None):
        P = self.poset
        if bottom is None and top is None:
            raise ValueError("Expected at least one element.")
        if bottom is None:
            bottom = P.bottom()
        if top is None and P.has_top():
            top = P.top()
        if top is None:
            PP = P.subposet(P.principal_order_filter(bottom))
        else:
            PP = P.subposet(P.closed_interval(bottom, top))
        return GradedPoset(poset=PP, R_label=self.R_label)

    def show(self):
        self.poset.show()

    @cached_method
    def Poincare_polynomial(self, abs_val=True):
        PR = PolynomialRing(QQ, 'Y')
        Y = PR.gens()[0]
        P = self.poset
        if abs_val:
            M = lambda x: abs(P.moebius_function(P.bottom(), x))
        else:
            M = lambda x: P.moebius_function(P.bottom(), x)
        return sum([M(x)*Y**P.rank_function()(x) for x in P])


def _combinatorial_eq_classes(GP:GradedPoset):
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
    return classes(GP.poset, GP.poset)