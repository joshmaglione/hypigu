#
#   Copyright 2020 Joshua Maglione
#
#   Distributed under MIT License
#

from functools import reduce as _reduce
from sage.misc.cachefunc import cached_method
from .Globals import __TIME as _time
from .Globals import __NCPUS as _N
import sage.parallel.decorate as _para


def _contract(M, rows):
    from sage.all import Matrix
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


# BUILD A WAY TO GET THE LABELS FROM THE CONTRACT MATRIX
def _get_labels(M, x, rows, L):
    from sage.all import VectorSpace, Set, Matrix

    # If non-central, then this is true iff not intersecting.
    not_e1 = lambda v: list(v)[1:] != [0]*(len(v)-1)

    # Determine new hyperplanes and group like rows together.
    V = VectorSpace(M.base_ring(), M.ncols())
    lines = []
    labels = []
    for r in range(M.nrows()):
        v = M[r]
        is_new = not_e1(v)
        i = 0
        while i < len(lines) and is_new:
            if v in V.subspace([lines[i]]):
                is_new = False
                labels[i] = labels[i].union(Set([r]))
            else:
                i += 1
        if is_new:
            lines.append(v)
            labels.append(Set([r]))

    # Adjust the labels because we are missing rows.
    fix_sets = lambda F: lambda S: Set(list(map(lambda s: F(s), S)))
    for k in rows:
        adjust = lambda i: i + (k <= i)*1
        labels = list(map(fix_sets(adjust), labels))
    labels = [Set(rows)] + labels

    # Adjust the row labels to hyperplane labels
    HL = L.hyperplane_labels
    A = L.hyperplane_arrangement
    HL_lab = lambda i: list(filter(lambda j: HL[j] == A[i], HL.keys()))[0]
    labels = list(map(fix_sets(HL_lab), labels))

    # Get the new hyperplanes
    FL = L.flat_labels
    new_hyp = [labels[k].union(labels[0]) for k in range(1, len(labels))]
    P = L.poset
    flat = lambda S: list(filter(lambda j: FL[j] == S, P.upper_covers(x)))[0]
    new_hyp = list(map(flat, new_hyp))

    def last_adj(S):
        l_hyp = labels[1:]
        T = S.difference(labels[0])
        new_T = Set([])
        for k in range(len(l_hyp)):
            if l_hyp[k].issubset(T):
                new_T = new_T.union(Set([new_hyp[k]]))
        return new_T

    return Matrix(lines), last_adj, {new_hyp[i]: i for i in range(len(new_hyp))}


def _parse_poset(P):
    global POS, atoms, labs, int_at
    from sage.all import Set
    import sage.parallel.decorate as para

    POS = P
    N = _N
    atoms = POS.upper_covers(POS.bottom())

    @para.parallel(N)
    def atom_set(k, shift):
        S = list(POS._elements[1+shift::k])
        m = lambda x: [x, Set(list(filter(lambda a: POS.le(a, x), atoms)))]
        return list(map(m, S))

    labs = list(atom_set([(N, k) for k in range(N)]))
    labs = [[P.bottom(), Set([])]] + _reduce(lambda x, y: x + y[1], labs, [])
    label_dict = {T[0]: T[1] for T in labs}

    return label_dict


# Can return the subarrangement A_x or restriction A^x simply based on the
# function F given. For subarrangement use 'lambda z: P.lower_covers(z)' and for
# restriction use 'lambda z: P.upper_covers(z)'.
def _subposet(P, x, F):
    from sage.all import Set
    elts = Set([])
    new_level = Set([x])
    while len(elts.union(new_level)) > len(elts):
        elts = elts.union(new_level)
        new_level = Set(_reduce(
            lambda x, y: x+y,
            map(F, new_level),
            []
        ))
    new_P = P.subposet(elts)
    return new_P


# Parallel function to build the intersection lattice.
# Moved to global to prevent accidentally carrying unnecessary data.
@_para.parallel(_N)
def build_next(A, S, HYP, LIN):
    from sage.all import exists, Set
    new_level = []
    new_hypcont = []
    if len(S) > 0:
        m = S[0]
    for i in S:
        T = LIN[i - m]
        for j in range(len(A)):
            # Skip the hyperplane already known to contain the intersection.
            if j not in HYP[i - m]:
                H = A[j]
                I = H._affine_subspace().intersection(T)
                # Check if the intersection is trivial.
                if I is not None:
                    if I == T:
                        # This case means that H cap T = T, so we should
                        # record that H contains T.
                        HYP[i - m] = HYP[i - m].union(Set([j]))
                    else:
                        # Check if we have this intersection already.
                        is_in, ind = exists(
                            range(len(new_level)),
                            lambda k: I == new_level[k]
                        )
                        if is_in:
                            # We have the intersection, so we update
                            # containment info accordingly.
                            new_hypcont[ind] = new_hypcont[ind].union(
                                Set([j]).union(HYP[i - m])
                            )
                        else:
                            # We do not have it, so we update everything.
                            new_level.append(I)
                            new_hypcont.append(HYP[i - m].union(Set([j])))
    return list(zip(new_level, new_hypcont))


# Default SageMath algorithm works well. However 'A.matroid()' seems to remove
# ordering, which we depend on, so care is needed.
def _lof_from_matroid(A=None, matroid=None):
    from sage.all import Set, Matrix, Matroid, Poset
    from functools import reduce
    if A is not None:
        rows = list(map(lambda H: H.coefficients()[1:], A.hyperplanes()))
        mat = Matrix(A.base_ring(), rows).transpose()
        M = Matroid(mat)
        n = len(A)
        lbl_map = lambda S: S
    else:
        M = matroid.simplify()
        n = len(M.groundset())
        grd_list = list(M.groundset())
        l_map = {x: grd_list.index(x) for x in grd_list}
        lbl_map = lambda S: frozenset([l_map[x] for x in S])

    L = M.lattice_of_flats()
    rank_r = lambda L, r: list(
        map(lbl_map, filter(lambda x: L.rank(x) == r, L._elements))
    )
    rank_1 = [frozenset([k]) for k in range(n)]
    ranks = reduce(
        lambda x, y: x + y,
        [rank_r(L, r) for r in range(2, L.rank() + 1)],
        [L.bottom()] + rank_1
    )
    P = Poset(
        L, element_labels={x: ranks.index(lbl_map(x)) for x in L._elements}
    )
    adj_set = lambda S: Set([x+1 for x in S])
    label_dict = {i: adj_set(ranks[i]) for i in range(len(L))}
    if A is not None:
        hyp_dict = {i: A[list(ranks[i])[0]] for i in range(1, n + 1)}
    else:
        hyp_dict = None
    return [P, label_dict, hyp_dict]


def _lof_from_affine_matroid(A):
    from sage.all import Poset, Set
    from functools import reduce
    A_coned = A.cone()
    hyps = list(map(lambda H: H.coefficients(), A_coned.hyperplanes()))
    extra = [0, 1] + [0]*(A.dimension())
    i = hyps.index(extra)
    P, L, H = _lof_from_matroid(A_coned)
    assert (H[i + 1]).coefficients() == extra
    new_elts = list(filter(lambda x: not P.le(i + 1, x), P))
    new_elts = reduce(
        lambda x, y: x + y,
        [list(filter(lambda x: P.rank(x) == r, new_elts)) for r in range(2, P.rank() + 1)],
        [0] + [j for j in range(1, len(A_coned) + 1) if j != i + 1]
    )
    new_names = {x: new_elts.index(x) for x in new_elts}
    adj_set = lambda S: Set([new_names[x] for x in S if x in new_elts])
    P_new = Poset(P.subposet(new_elts), element_labels=new_names)
    L_new = {new_names[x]: adj_set(L[x]) for x in new_elts}

    def inv_H(h):
        cut = lambda x: x.coefficients()[1:]
        pair = list(filter(lambda x: cut(x[1]) == h, H.items()))
        return pair[0][0]
    H_new = {new_names[inv_H(h.coefficients())]: h for h in A}
    return [P_new, L_new, H_new]


class LatticeOfFlats():

    def __init__(self, A=None, poset=None, flat_labels=None,
    hyperplane_labels=None, lazy=False, matroid=None,
    nature_hyperplane_label=True):
        self.hyperplane_arrangement = A
        self.poset = poset
        self.flat_labels = flat_labels
        self.hyperplane_labels = hyperplane_labels
        if poset is not None:
            assert poset.has_bottom(), "Expected a unique minimal element in poset."
            assert poset.is_graded(), "Expected a graded poset."
            self.poset = poset
        else:
            if not lazy:
                if A is not None:
                    if A.is_central():
                        P, FL, HL = _lof_from_matroid(A)
                    else:
                        P, FL, HL = _lof_from_affine_matroid(A)
                else:
                    P, FL, HL = _lof_from_matroid(A=None, matroid=matroid)
                self.poset = P
                self.flat_labels = FL
                self.hyperplane_labels = HL
        if self.flat_labels is None and not lazy:
            self.flat_labels = _parse_poset(poset)
        if self.hyperplane_arrangement is not None and self.hyperplane_labels is None and nature_hyperplane_label:
            self.hyperplane_labels = {i + 1: A[i] for i in range(len(A))}

    def __repr__(self):
        if self.hyperplane_arrangement:
            return "The lattice of flats of:\n{0}\ngiven by:\n{1}".format(self.hyperplane_arrangement, self.poset)
        else:
            return "The lattice of flats of some matroid given by:\n{0}".format(self.poset)

    def _save(self, file, var_name='L'):
        from sage.all import Matrix
        HH = self.hyperplane_arrangement.parent()
        A = Matrix(map(lambda H: H.coefficients(), self.hyperplane_arrangement.hyperplanes())).rows()
        CR = tuple(map(lambda T: tuple(T), self.poset.cover_relations()))
        FL = self.flat_labels
        FL_tup = tuple([tuple([x, list(FL[x])]) for x in FL.keys()])
        del FL
        dict_builder = "FL = {x[0]: Set(x[1]) for x in FL_tup}\n"
        with open(file, "w") as F:
            F.write("from sage.all import HyperplaneArrangements, QQ, Poset, Set\n")
            F.write("import hypigu as hi\n")
            F.write("H = HyperplaneArrangements(QQ, {0})\n".format(HH.variable_names()))
            del HH
            F.write("A = H({0})\n".format(A).replace("), ", "),\n"))
            del A
            F.write("CR = {0}\n".format(CR).replace("), ", "),\n"))
            del CR
            F.write("P = Poset([range({0}), CR], cover_relations=True)\n".format(len(self.poset._elements)))
            F.write("FL_tup = {0}\n".format(FL_tup).replace("), ", "),\n"))
            F.write(dict_builder)
            F.write("del FL_tup\n")
            F.write("{0} = hi.LatticeOfFlats(A, poset=P, flat_labels=FL)\n".format(var_name))
            F.write("del H, A, CR, P, FL\n")
            F.write("print('Loaded a lattice of flats. Variable name: {0}')".format(var_name))

    def atoms(self):
        return self.poset.upper_covers(self.poset.bottom())

    def labels_of_flats(self):
        elt_tup = lambda x: tuple([x, self.flat_labels[x]])
        return list(map(elt_tup, self.poset._elements))

    def labels_of_hyperplanes(self):
        P = self.poset
        elt_tup = lambda x: tuple([x, self.hyperplane_labels[x]])
        return list(map(elt_tup, P.upper_covers(P.bottom())))

    def proper_part_poset(self):
        P = self.poset
        elts = list(P._elements)
        if P.has_top():
            elts.remove(P.top())
        elts.remove(P.bottom())
        return P.subposet(elts)

    def show(self):
        self.poset.show()

    def subarrangement(self, x):
        P = self.poset
        if type(x) != set:
            assert x in P, "Expected element to be in poset."
            new_P = _subposet(P, x, lambda z: P.lower_covers(z))
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
            return LatticeOfFlats(new_A, poset=new_P, flat_labels=new_FL, hyperplane_labels=new_HL)
        else:
            L = self.flat_labels
            X = list(filter(lambda y: L[y] == x, P._elements))
            try:
                return self.subarrangement(X[0])
            except IndexError:
                raise ValueError("No element labeled by:\n{0}".format(x))

    def restriction(self, x):
        from sage.all import Matrix, HyperplaneArrangements
        P = self.poset
        if type(x) != set:
            assert x in P, "Expected element to be in poset."
            new_P = _subposet(P, x, lambda z: P.upper_covers(z))
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
            return LatticeOfFlats(new_A, poset=new_P, flat_labels=new_FL, hyperplane_labels=new_HL)
        else:
            L = self.flat_labels
            X = list(filter(lambda y: L[y] == x, P._elements))
            try:
                return self.restriction(X[0])
            except IndexError:
                raise ValueError("No element labeled by:\n{0}".format(x))

    def deletion(self, H):
        from sage.all import Set

        P = self.poset
        L = self.flat_labels

        if type(H) == set:
            L = self.flat_labels
            X = list(filter(lambda y: L[y] == H, P._elements))
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
            return bool(len(S) == m and H not in S)
        new_top = list(filter(check, coatoms))

        if len(new_top) == 1:
            new_P = _subposet(P, new_top[0], lambda z: P.lower_covers(z))
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

        return LatticeOfFlats(new_HPA, poset=new_P, flat_labels=new_FL, hyperplane_labels=new_HL)

    def _lazy_restriction(self, H):
        HPA = self.hyperplane_arrangement
        assert HPA is not None, "Needs underlying hyperplane arrangement."
        A = HPA.restriction(HPA[H - 1])
        return LatticeOfFlats(A, lazy=True)

    def _lazy_deletion(self, H):
        HPA = self.hyperplane_arrangement
        assert HPA is not None, "Needs underlying hyperplane arrangement."
        return LatticeOfFlats(HPA.parent()(HPA[:H-1] + HPA[H:]), lazy=True)

    @cached_method
    def Poincare_polynomial(self, abs_val=False):
        from sage.all import QQ, PolynomialRing
        if abs_val not in {True, False}:
            raise ValueError("abs must be either True or False.")
        PR = PolynomialRing(QQ, 'Y')
        Y = PR.gens()[0]
        if self.poset is not None:
            P = self.poset
            atoms = self.atoms()
            if P.rank() == 0:
                return PR(1)
            if P.rank() == 1:
                return PR(1 + len(atoms)*Y)
        else:
            # Lazy
            A = self.hyperplane_arrangement
            assert A is not None, "Expected either a poset or hyperplane arrangement."
            if A.rank() == 0:
                return PR(1)
            if A.rank() == 1:
                return PR(1 + len(A)*Y)
        if self.hyperplane_arrangement is not None:
            try:  # Some hyperplane arrangements are bugged in SageMath.
                D = self._lazy_deletion(1)
                R = self._lazy_restriction(1)
                return PR(D.Poincare_polynomial() + Y*R.Poincare_polynomial())
            except Exception:
                pass
        if abs_val:
            comb_elts = self._combinatorial_eq_elts()
            rk = P.rank_function()
            bot = P.bottom()
            pi = PR(1) 
            if P.has_top():
                top = P.top()
                pi += abs(P.moebius_function(bot, top))*Y**(rk(top))
            for t in comb_elts:
                pi += abs(P.moebius_function(bot, t[0]))*t[1]*Y**(rk(t[0]))
            return pi 
        else:
            chi = P.characteristic_polynomial()
            q = chi.variables()[0]
            d = chi.degree(q)
            return PR((-Y)**d*chi.subs({q: -Y**-1}))

    @cached_method
    def _combinatorial_eq_elts(self):
        global POS, P_elts
        import sage.parallel.decorate as para

        N = _N
        POS = self.poset
        P_elts = self.proper_part_poset()._elements

        @para.parallel(N)
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
        prelim_elts = _reduce(lambda x, y: x + y[1], prelim_elts, [])

        # Test further to minimize the size.
        equiv_elts = []
        while len(prelim_elts) > 0:
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

