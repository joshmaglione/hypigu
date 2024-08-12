#
#   Copyright 2020--2024 Joshua Maglione
#
#   Distributed under MIT License
#
from sage.all import Matrix, VectorSpace, Set, exists, HyperplaneArrangements, Poset, Matroid, QQ, PolynomialRing, DiGraph

from functools import reduce
from sage.misc.cachefunc import cached_method
from .globals import _TIME as _time
from .globals import _NCPUS as _N
import sage.parallel.decorate as _para


def _contract(M, rows):
    K = M.base_ring()
    Q = [M[k] for k in rows] + [M[k] for k in range(M.nrows()) if k not in rows]
    Q = Matrix(K, Q).transpose()
    top = Q[0]
    E = Q[1:, :].echelon_form()
    rel_E = E[:, :len(rows)]
    out_E = E[:, len(rows):]
    rm = rel_E.pivot_rows()
    # In case we have a non-central arrangement
    for r in rm:
        v = E[r]
        i = list(v).index(1)
        top = top - top[i]*v
    A = [tuple(list(top)[len(rows):])] + [out_E[k] for k in range(E.nrows())
                                          if k not in rm]
    M_out = Matrix(K, A).transpose()
    return M_out

# Can return the subarrangement A_x or restriction A^x simply based on the
# function F given. For subarrangement use 'P.lower_covers' and for
# restriction use 'P.upper_covers'.
def _subposet(P, x, F):
    elts = Set()
    new_level = Set([x])
    while len(elts.union(new_level)) > len(elts):
        elts = elts.union(new_level)
        new_level = Set(reduce(
            lambda x, y: x+y,
            map(F, new_level),
            []
        ))
    return P.subposet(elts)


class GradedPoset():

    def __init__(
            self,
            arrangement=None,
            matroid=None,
            poset=None,
        ):
        self.arrangement = arrangement
        self.matroid = matroid
        self.poset = poset
        if arrangement is None and matroid is None and poset is None:
            raise ValueError("Expected either an arrangement, matroid, or poset.")
        if poset is not None:
            if not poset.has_bottom():
                raise ValueError("Expected a unique minimal element in poset.")
            if not poset.is_graded(): 
                raise ValueError("Expected a graded poset.")
        if matroid is not None:
            self.poset = matroid.lattice_of_flats()
        if arrangement is not None:
            M = Matrix([H.normal() for H in arrangement])
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
        start = f"A graded poset with {len(self.poset)} elements"
        if self.arrangement:
            return f"{start} built from the intersection poset of\n{self.arrangement}"
        if self._matroid:
            return f"{start} built from the lattice of flats of\n{self.matroid}"
        return start

    def atoms(self):
        return self.poset.upper_covers(self.poset.bottom())

    def show(self):
        self.poset.show()

    def subarrangement(self, x):
        P = self.poset
        if not isinstance(x, set):
            assert x in P, "Expected element to be in poset."
            new_P = _subposet(P, x, P.lower_covers)
            new_A = None
            new_FL = None
            new_HL = None
            if self.hyperplane_arrangement and self.hyperplane_labels:
                A = self.hyperplane_arrangement
                HL = self.hyperplane_labels
                atoms = new_P.upper_covers(new_P.bottom())
                keep = list(map(lambda k: HL[k], atoms))
                new_A = A.parent()(keep)
                new_HL = {a: HL[a] for a in atoms}
            if self.flat_labels:
                FL = self.flat_labels
                new_FL = {x: FL[x] for x in new_P._elements}
            return GradedPoset(new_A, poset=new_P, flat_labels=new_FL, hyperplane_labels=new_HL)
        else:
            L = self.flat_labels
            X = list(filter(lambda y: L[y] == x, P._elements))
            try:
                return self.subarrangement(X[0])
            except IndexError:
                raise ValueError("No element labeled by:\n{0}".format(x))

    def restriction(self, x):
        P = self.poset
        if not isinstance(x, set):
            assert x in P, "Expected element to be in poset."
            new_P = _subposet(P, x, P.upper_covers)
            new_A = None
            new_HL = None
            if self.hyperplane_arrangement:
                A = self.hyperplane_arrangement
                hyp_coeffs = map(lambda H: H.coefficients(), A.hyperplanes())
                M = Matrix(A.base_ring(), list(hyp_coeffs))
                rows = sorted(list(map(
                    lambda H: list(A).index(self.hyperplane_labels[H]),
                    self.flat_labels[x]
                )))
                new_M = _contract(M, rows)
                new_M, lab_func, hyp_dict = _get_labels(new_M, x, rows, self)
                HH = HyperplaneArrangements(
                    A.base_ring(),
                    A.parent().variable_names()[:new_M.ncols()-1]
                )
                new_A = HH(new_M)
                FL = self.flat_labels
                new_FL = {x: lab_func(FL[x]) for x in new_P._elements}
                new_HL = {a: new_A[hyp_dict[a]] for a in new_P.upper_covers(new_P.bottom())}
            else:
                FL = self.flat_labels
                new_FL = {y: FL[y].difference(FL[x]) for y in new_P._elements}
            return GradedPoset(new_A, poset=new_P, flat_labels=new_FL, hyperplane_labels=new_HL)
        else:
            L = self.flat_labels
            X = list(filter(lambda y: L[y] == x, P._elements))
            try:
                return self.restriction(X[0])
            except IndexError:
                raise ValueError("No element labeled by:\n{0}".format(x))

    def deletion(self, H):
        P = self.poset
        L = self.flat_labels
        if isinstance(H, set):
            L = self.flat_labels
            X = [y for y in P._elements if L[y] == H]
            try:
                return self.deletion(X[0])
            except IndexError:
                raise ValueError("No element labeled by:\n{0}".format(H))

        assert P.rank_function()(H) == 1, "Expected an atom."

        if P.has_top():
            coatoms = P.lower_covers(P.top())
        else:
            # not really coatoms... but whatever
            coatoms = P.maximal_elements()

        m = len(self.atoms()) - 1

        def check(C):
            S = L[C]
            return len(S) == m and H not in S

        new_top = list(filter(check, coatoms))

        if len(new_top) == 1:
            new_P = _subposet(P, new_top[0], P.lower_covers)
            new_FL = {y: L[y] for y in new_P._elements}
        else:
            def good_flats(F):
                S = L[F]
                if H in S:
                    U = S.difference(Set([H]))
                    return (U != 0) and (U not in L.values())
                else:
                    return True
            flats = list(filter(good_flats, P._elements))
            new_P = P.subposet(flats)
            new_FL = {y: L[y].difference(Set([H])) for y in flats}

        if self.hyperplane_arrangement:
            HPA = self.hyperplane_arrangement
            HL = self.hyperplane_labels
            A = list(HPA)
            A.remove(HPA[H])
            new_HPA = HPA.parent()(A)
            new_HL = {x: HL[x] for x in new_P.upper_covers(new_P.bottom())}
        else:
            new_HPA = None
            new_HL = None

        return GradedPoset(new_HPA, poset=new_P, flat_labels=new_FL, hyperplane_labels=new_HL)

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

    @cached_method
    def _combinatorial_eq_elts(self):
        global POS, P_elts

        N = _N
        POS = self.poset
        P_elts = self.proper_part_poset()._elements

        @_para.parallel(N)
        def match_elts(k, shift):
            all_elts = P_elts[shift::k]
            eq_elts = []
            counts = []
            down = []
            restrict = []
            while len(all_elts) > 0:
                x = all_elts[0]
                dow_x = self.subarrangement(x)
                res_x = self.restriction(x)
                match = False
                i = 0
                while not match and i < len(eq_elts):
                    if dow_x.poset.is_isomorphic(down[i].poset) and res_x.poset.is_isomorphic(restrict[i].poset):
                        match = True
                    else:
                        i += 1
                if match:
                    counts[i] += 1
                else:
                    eq_elts.append(x)
                    counts.append(1)
                    down.append(dow_x)
                    restrict.append(res_x)
                all_elts = all_elts[1:]
            return list(zip(eq_elts, counts, down, restrict))

        # Get the preliminary set of inequivalent elements
        prelim_elts = list(match_elts([(N, k) for k in range(N)]))
        prelim_elts = reduce(lambda x, y: x + y[1], prelim_elts, [])

        # Test further to minimize the size.
        equiv_elts = []
        while prelim_elts:
            x = prelim_elts[0]
            match = False
            i = 0
            while not match and i < len(equiv_elts):
                if x[2].poset.is_isomorphic(equiv_elts[i][2].poset) and x[3].poset.is_isomorphic(equiv_elts[i][3].poset):
                    match = True
                else:
                    i += 1
            if match:
                equiv_elts[i][1] += x[1]
            else:
                equiv_elts.append(list(x))
            prelim_elts = prelim_elts[1:]
        return equiv_elts
