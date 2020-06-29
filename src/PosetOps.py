#
#   Copyright 2020 Joshua Maglione 
#
#   Distributed under MIT License
#

def _my_intersection_poset(A):
    from sage.geometry.hyperplane_arrangement.affine_subspace import AffineSubspace
    from sage.all import flatten, QQ, VectorSpace
    K = A.base_ring()
    from sage.geometry.hyperplane_arrangement.affine_subspace import AffineSubspace
    whole_space = AffineSubspace(0, VectorSpace(K, A.dimension()))
    L = [[whole_space]]
    label = [[[]]]
    active = True
    codim = 0
    while active:
        active = False
        new_level = []
        new_label = []
        for i in range(len(L[codim])):
            T = L[codim][i]
            for j in range(len(A)):
                H = A[j]
                I = H._affine_subspace().intersection(T)
                if I is not None and I != T and I not in new_level:
                    new_level.append(I)
                    new_label.append(label[codim][i] + [j])
                    active = True
        if active:
            L.append(new_level)
            label.append(new_label)
        codim += 1
    from sage.misc.flatten import flatten
    L = flatten(L)
    label = reduce(lambda x, y: x + y, label, [])
    t = {}
    for i in range(len(L)):
        t[i] = L[i]
    cmp_fn = lambda p, q: t[q] < t[p]
    list_str = lambda L : reduce(lambda x, y: x + ' ' + str(y), L[1:], str(L[0]))
    elt_labels = [''] + list(map(list_str, label[1:]))
    from sage.combinat.posets.posets import Poset
    return Poset(
        (t, cmp_fn), 
        element_labels=elt_labels
    )

def CharacteristicFunction(A):
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
    P = _my_intersection_poset(A)
    def upper_subposet(X):
        S = set(P.upper_covers(X))
        checked = set()
        while len(S - checked) > 0:
            x = (S - checked).pop()
            S = S.union(set(P.upper_covers(x)))
            checked.add(x)
        return P.subposet([X] + list(S))
    def char_func(X):
        PX = upper_subposet(X)
        moebius = map(lambda Y: PX.moebius_function(X, Y), PX._elements)
        bot = P.bottom()
        my_rank = lambda Y: P.subposet(P.closed_interval(bot, Y)).height() - 1
        dims = map(lambda Y: A.dimension() - my_rank(Y), PX._elements)
        t = var('x')
        return reduce(lambda x, y: x + y[0]*t**y[1], zip(moebius, dims), 0)
    return char_func, P