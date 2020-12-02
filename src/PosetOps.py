#
#   Copyright 2020 Joshua Maglione 
#
#   Distributed under MIT License
#

from functools import reduce as _reduce

# Returns the proper part of a poset as a subposet.
def _proper_part(P):
    central = P.has_top()
    if central: 
        prop = filter(lambda X: X != P.top() and X != '', P._elements)
    else: 
        prop = filter(lambda X: X != '', P._elements)
    return P.subposet(prop)


# Can return the subarrangement A_x or restriction A^x simply based on the
# function F given. For deletion use 'lambda z: P.lower_covers(z)' and for
# restriction use 'lambda z: P.upper_covers(z)'.
def _subposet(P, x, F):
    from sage.all import Set, Poset, DiGraph
    elts = Set([])
    if type(x) == list:
        new_level = Set(x)
    else:
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

# Returns the characteristic polynomial of P -- using the fact that it comes
# from a hyperplane arrangement. 
# Copied from the Hyperplane Arrangement code in Sage 9.2.
def char_poly(P, dim=0):
    from sage.all import QQ
    from sage.rings.polynomial.polynomial_ring import polygen

    if dim == 0:
        dim = P.rank() 
    atoms = P.upper_covers(P.bottom())
    X = polygen(QQ, 'X')
    if P.rank() == 2 and len(atoms) == 2: 
        return X**dim - 2*X**(dim - 1) + X**(dim - 2)
    if P.rank() == 1:
        return X**(dim - 1)*(X - len(atoms))
    H = atoms[0]
    R = _subposet(P, H, lambda z: P.upper_covers(z))
    coatoms = P.lower_covers(P.top())
    sub_coat = list(filter(lambda x: not P.le(H, x), coatoms))
    D = P.subposet(list(_subposet(P, sub_coat, lambda z: P.lower_covers(z))._elements) + [P.top()])
    charpoly_R = char_poly(R, dim=dim - 1)
    charpoly_D = char_poly(D, dim=dim)
    return charpoly_D - charpoly_R

# Two elements x, y of P are equivalent if A_x = A_y and A^x = A^y. We need only
# a representative of each equivalence class and some other data. We return a
# list of triples: a rep. elt, size of equiv. class, and the poset from A_x.
def _equiv_elts(P):
    global POS, P_elts
    from os import cpu_count
    import sage.parallel.decorate as para

    N = cpu_count()
    POS = P
    P_elts = POS._elements[1:-1] # Proper part

    @para.parallel(N)
    def match_elts(k, shift):
        all_elts = P_elts[shift::k]
        eq_elts = []
        counts = []
        deletion = []
        restrict = []
        while len(all_elts) > 0:
            x = all_elts[0]
            del_x = _subposet(POS, x, lambda z: POS.lower_covers(z))
            res_x = _subposet(POS, x, lambda z: POS.upper_covers(z))
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
        return list(zip(eq_elts, counts, deletion, restrict))

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
            if x[2].is_isomorphic(equiv_elts[i][2]) and x[3].is_isomorphic(equiv_elts[i][3]):
                match = True
            else:
                i += 1
        if match:
            equiv_elts[i][1] += x[1]
        else:
            equiv_elts.append(list(x))
        prelim_elts = prelim_elts[1:]
    return equiv_elts

# The deletion method is coded well enough to allow for this kind of recursion. 
# We need the new labels for a hyperplane though.
def _deletion(A, X, P, poset=True, OG=None):
    if OG:
        Ar = OG
    else:
        Ar = A
    H_to_str = lambda H: str(list(Ar).index(H))
    complement = filter(lambda H: not P.le(H_to_str(H), X), A)
    B = _reduce(lambda x, y: x.deletion(y), complement, A)
    if not poset:
        return B
    Q = _subposet(P, X, lambda z: P.lower_covers(z))
    return B, Q


def _Coxeter_poset_data():
    # Bell numbers: A000110
    def A_poset(n):
        from sage.all import binomial
        S = [1, 1, 2, 5, 15, 52, 203]
        while len(S) <= n+1:
            m = len(S) - 1
            S.append(_reduce(
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
        from sage.all import add
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



# We expand on the function in sage, optimizing a little bit. This makes little
# difference in small ranks but noticeable difference in larger ranks. This is
# still quite slow. 
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
    list_str = lambda L : _reduce(lambda x, y: x+' '+str(y), L[1:], str(L[0]))
    elt_labels = [''] + list(map(list_str, L[1:]))
    
    return Poset((t, cmp_fn), element_labels=elt_labels)

# Compute the Poincare polynomial of either the chain F or the upper ideal of P
# at restrict. 
def PoincarePolynomial(P, F=None, restrict=None):
    from .Globals import __DEFAULT_p as p
    from sage.all import var, ZZ

    # Setup the data
    if restrict != None:
        F = [restrict]
    if F == None:
        assert P.has_bottom()
        F = [P.bottom()]
    if F[0] != P.bottom():
        F = [P.bottom()] + F 

    P_up = _subposet(P, F[-1], lambda z: P.upper_covers(z))
    chi = P_up.characteristic_polynomial() 
    d = chi.degree()
    p = var(p)
    pi = ((-p)**d*chi(q=-p**(-1))).expand().simplify()
    if F[-1] == P.bottom() or restrict != None:
        return pi
    P_down = _subposet(P, F[-1], lambda z: P.lower_covers(z))
    return pi*PoincarePolynomial(P_down, F=F[:-1])
