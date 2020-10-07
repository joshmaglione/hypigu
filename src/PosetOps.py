#
#   Copyright 2020 Joshua Maglione 
#
#   Distributed under MIT License
#

# Returns the proper part of a poset as a subposet.
def _proper_part(P):
    central = P.has_top()
    if central: 
        prop = filter(lambda X: X != P.top() and X != '', P._elements)
    else: 
        prop = filter(lambda X: X != '', P._elements)
    return P.subposet(prop)

# Labels the n coatoms of P with {0, ..., n-1}, and updates everything above.
def _update_labels(P):
    from sage.all import Poset, DiGraph, Set, exists
    hyps = P.upper_covers('')
    new_labels = [str(k) for k in range(len(hyps))]
    hyp_dict = {hyps[i] : new_labels[i] for i in range(len(hyps))}

    def update_elt(x):
        if x == '':
            return ''
        if not ' ' in x:
            return hyp_dict[x]
        hyps_below = filter(lambda H: P.le(H, x), hyps)
        new_x_labels = map(lambda H: hyp_dict[H], hyps_below)
        return reduce(lambda u, v: u + v + ' ', new_x_labels, '')

    digraph_dat = map(
        lambda x: [update_elt(x[0]), update_elt(x[1])], 
        P.cover_relations()
    )
    G = DiGraph(digraph_dat)
    assert Poset(G).is_isomorphic(P)
    return Poset(G)

# Can return the deletion A_x or restriction A^x simply based on the function F
# given. For deletion use 'lambda z: P.lower_covers(z)' and for restriction use
# 'lambda z: P.upper_covers(z)'.
def _subposet(P, x, F, update=False):
    from sage.all import Set, Poset, DiGraph
    elts = Set([])
    new_level = Set([x])
    while len(elts.union(new_level)) > len(elts):
        elts = elts.union(new_level)
        new_level = Set(reduce(
            lambda x, y: x+y, 
            map(F, new_level), 
            []
        ))
    new_P = Poset(DiGraph(P.subposet(elts).cover_relations()))
    assert new_P.is_isomorphic(P.subposet(elts))
    if update:
        return _update_labels(new_P)
    return new_P

# Two elements x, y of P are equivalent if A_x = A_y and A^x = A^y. We need only
# a representative of each equivalence class and some other data. We return a
# list of triples: a rep. elt, size of equiv. class, and the poset from A_x.
def _equiv_elts(P):
    P_prop = _proper_part(P)
    all_elts = P_prop._elements
    eq_elts = []
    counts = []
    deletion = []
    restrict = []
    while len(all_elts) > 0:
        x = all_elts[0]
        del_x = _subposet(P, x, lambda z: P.lower_covers(z), update=True)
        res_x = _subposet(P, x, lambda z: P.upper_covers(z))
        match = False
        i = 0
        while not match and i < len(eq_elts):
            if del_x.is_isomorphic(deletion[i]) and res_x.is_isomorphic(restrict[i]):
                match = True
            else:
                i += 1
        if match:
            counts[i] += 1
        else:
            eq_elts.append(x)
            counts.append(1)
            deletion.append(del_x)
            restrict.append(res_x)
        all_elts = all_elts[1:]
    return zip(eq_elts, counts, deletion) 


# The deletion method is coded well enough to allow for this kind of recursion. 
# We need the new labels for a hyperplane though.
def _deletion(A, X, P):
    H_to_str = lambda H: str(list(A).index(H))
    complement = filter(lambda H: not P.le(H_to_str(H), X), A)
    B = reduce(lambda x, y: x.deletion(y), complement, A)
    # A dictionary to convert old labels to new labels.
    def new_labels(k):
        for i in range(len(B)):
            if B[i] == A[k]:
                return i
        return -1
    return B, {str(k): str(new_labels(k)) for k in range(len(A))}


# We expand on the function in sage, optimizing a little bit. This makes little
# difference in small ranks but noticeable different in larger ranks. 
def IntersectionPoset(A):
    from sage.geometry.hyperplane_arrangement.affine_subspace import AffineSubspace
    from sage.all import exists, flatten, Set, QQ, VectorSpace, Poset
    K = A.base_ring()
    whole_space = AffineSubspace(0, VectorSpace(K, A.dimension()))
    # L is the ranked list of affine subspaces in L(A).
    L = [[whole_space], list(map(lambda H: H._affine_subspace(), A))]
    # hyp_cont is the ranked list describing which hyperplanes contain the
    # corresponding intersection. 
    hyp_cont = [[Set([])], [Set([k]) for k in range(len(A))]]
    active = True
    codim = 1
    while active:
        active = False
        new_level = []
        # new_label = []
        new_hypcont = []
        for i in range(len(L[codim])):
            T = L[codim][i]
            for j in range(len(A)):
                # Skip the hyperplane already known to contain the intersection.
                if not j in hyp_cont[codim][i]: 
                    H = A[j]
                    I = H._affine_subspace().intersection(T)
                    # Check if the intersection is trivial.
                    if I is not None:
                        if I == T: 
                            # This case means that H cap T = T, so we should
                            # record that H contains T.
                            hyp_cont[codim][i] = hyp_cont[codim][i].union(
                                Set([j])
                            )
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
                                    Set([j]).union(hyp_cont[codim][i])
                                )
                            else:
                                # We do not have it, so we update everything.
                                new_level.append(I)
                                new_hypcont.append(
                                    hyp_cont[codim][i].union(Set([j]))
                                )
                                active = True
        if active:
            L.append(new_level)
            hyp_cont.append(new_hypcont)
        codim += 1
    
    L = flatten(hyp_cont)
    t = {}
    for i in range(len(L)):
        t[i] = L[i]
    cmp_fn = lambda p, q: t[p].issubset(t[q])
    list_str = lambda L : reduce(lambda x, y: x+' '+str(y), L[1:], str(L[0]))
    elt_labels = [''] + list(map(list_str, L[1:]))
    
    return Poset((t, cmp_fn), element_labels=elt_labels)


def CharacteristicFunction(A, poset=None):
    r"""
    Given a hyperplane arrangement A, return a function chi on the elements of 
    the intersection poset of A such that chi(X) returns a polynomial in ``x``
    equal to the characteristic polynomial of A^X. Here, the elements of the
    poset are assumed to be ``range(len(A.intersection_poset()))``. 

    INPUT:

    - ``A`` -- hyperplane arrangement. 

    OUTPUT: a function on range(len(P)), where P is the intersection poset of H.

    EXAMPLES:

    This example illustrates how to build a Coxeter arrangement ::

        sage: H = CoxeterArrangement("B", 4)
        sage: H
        Arrangement of 16 hyperplanes of dimension 4 and rank 4

    We can also combine the rank into the string as follows ::

        sage: H = IA.CoxeterArrangement("B4")
        sage: H
        Arrangement of 16 hyperplanes of dimension 4 and rank 4

    """
    from sage.all import var
    if not poset:
        P = IntersectionPoset(A)
    else:
        P = poset
    def upper_subposet(X):
        return _subposet(P, X, lambda z: P.upper_covers(z))
    def char_func(X):
        from Globals import __DEFAULT_p as p
        p = var(p)
        PX = upper_subposet(X)
        moebius = map(lambda Y: PX.moebius_function(X, Y), PX._elements)
        bot = P.bottom()
        # my_rank = lambda Y: P.subposet(P.closed_interval(bot, Y)).height() - 1
        dims = map(lambda Y: A.dimension() - P.rank(Y), PX._elements) # codim
        return reduce(lambda x, y: x + y[0]*p**y[1], zip(moebius, dims), 0)
    return char_func, P


def PoincarePolynomial(A, F=[''], poset=None):
    from Globals import __DEFAULT_p as p
    from sage.all import var, ZZ
    p = var(p)
    char_func, P = CharacteristicFunction(A, poset=poset)
    chi = char_func(F[-1])
    d = chi.degree(p)
    pi = ((-p)**d*chi.subs({p:-p**-1})).factor().simplify()
    if F[-1] == '':
        return pi
    # X = map(lambda i: ZZ(i), F[-1].split(' '))
    X = F[-1]
    B, new_labels = _deletion(A, X, P)
    def update_str(Y):
        if Y == '':
            return Y 
        L = Y.split(' ')
        L_upd = map(lambda s: new_labels[s], L)
        return reduce(lambda x, y: x + ' ' + y, L_upd[1:], L_upd[0])
    G = map(update_str, F[:-1])
    return pi*PoincarePolynomial(B, G)




def _Coxeter_poset_data():
    # Bell numbers: A000110
    def A_poset(n):
        from sage.all import binomial
        S = [1, 1, 2, 5, 15, 52, 203]
        while len(S) <= n+1:
            m = len(S) - 1
            S.append(reduce(
                lambda x, y: x + y[0]*binomial(m, y[1]), zip(S, range(m+1)), 0
            ))
        return S[n+1]
    def S(n, k, m):
        if k > n or k < 0 : 
            return 0
        if n == 0 and k == 0: 
            return 1
        return S(n-1, k-1, m) + (m*(k+1)-1)*S(n-1, k, m)
    def A007405(n): 
        return add(S(n, k, 2) for k in (0..n)) # Peter Luschny, May 20 2013
    # D analog of Bell numbers: A039764
    Dlist = [1, 1, 4, 15, 72, 403, 2546, 17867, 137528, 1149079, 10335766, 99425087, 1017259964, 11018905667, 125860969266, 1510764243699, 18999827156304, 249687992188015, 3420706820299374, 48751337014396167]
    table = {
        'A': {
            'hyperplanes': lambda n: n*(n+1) // 2,
            'poset': A_poset
        },
        'B': {
            'hyperplanes': lambda n: n**2,
            'poset': A007405
        },
        'D': {
            'hyperplanes': lambda n: n**2 - n,
            'poset': lambda n: Dlist[n]
        }
    }
    return table

def _possibly_Coxeter(P):
    r = P.rank()
    hypers = list(filter(lambda x: P.covers(P.bottom(), x), P))
    m = len(hypers)
    CPD = _Coxeter_poset_data()
    for name in ['A', 'B', 'D']:
        if CPD[name]['hyperplanes'](r) == m:
            if CPD[name]['poset'](r) == len(P):
                return [True, name]
    return [False, None]


