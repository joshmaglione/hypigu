#
#   Copyright 2020 Joshua Maglione 
#
#   Distributed under MIT License
#

# A file to profile some of the sage functions against my own.

# A and B are affine subspaces
def _Aff_intersection(A, B):
    import sympy as sp
    from sage.all import SR, Matrix, QQ, VectorSpace
    from sage.geometry.hyperplane_arrangement.affine_subspace import AffineSubspace
    # Setup
    p_A = A.point()
    p_B = B.point()
    M_A = A.linear_part().matrix()
    M_B = B.linear_part().matrix()
    d_A = M_A.nrows()
    d_B = M_B.nrows()
    n = M_A.ncols()
    M = sp.Matrix(M_A.rows() + (-M_B).rows()).transpose()
    b = list(p_B - p_A)
    col_b = sp.Matrix(n, int(1), b)
    symbs = sp.symbols('x0:%s' % (d_A+d_B))
    V = VectorSpace(QQ, n)

    # Solves
    Mat = M.rref()[0]
    N = Mat.nullspace()
    v = sp.linsolve(tuple([Mat, col_b]), symbs)
    if len(v) == 0:
        return None
    v = next(iter(v))

    # Massage the data
    get_vect = lambda y: Matrix(1, d_A, 
        map(lambda x: QQ(x._sage_()), list(y)[:d_A])
    )
    mult_M_A = lambda x: (x*M_A).list()
    eval_func = lambda f: SR(f).subs({SR(s) : 0 for s in symbs})
    null_basis = map(get_vect, N)
    pt = V(mult_M_A(get_vect(map(eval_func, v)))) + p_A
    basis = list(map(mult_M_A, null_basis))

    return AffineSubspace(pt, V.subspace(basis))

def _Aff_eq(A, B):
    if A.dimension() != B.dimension():
        return False
    else:
        I = _Aff_intersection(A, B)
        if I: 
            return bool(I.dimension() == A.dimension())
        return False

def profile(N):
    from sage.geometry.hyperplane_arrangement.affine_subspace import AffineSubspace
    from sage.all import randint, random_matrix, random_vector, QQ, VectorSpace, ZZ
    import time
    dims = [randint(3, 20) for _ in range(N)]
    sub_dim1 = map(lambda x: randint(1, x), dims)
    sub_dim2 = map(lambda x: randint(1, x), dims)
    vsps = map(lambda d: VectorSpace(QQ, d), dims)
    dim_to_pt = lambda d: random_vector(QQ, d)
    dims_to_sub = lambda d: d[0].subspace(random_matrix(QQ, d[2], d[1]).rows())
    pts1 = map(dim_to_pt, dims)
    pts2 = map(dim_to_pt, dims)
    sub1 = map(dims_to_sub, zip(vsps, dims, sub_dim1))
    sub2 = map(dims_to_sub, zip(vsps, dims, sub_dim2))
    Aff1 = list(map(lambda x: AffineSubspace(x[0], x[1]), zip(pts1, sub1)))
    Aff2 = list(map(lambda x: AffineSubspace(x[0], x[1]), zip(pts2, sub2)))

    # Equality test:
    my_eq = lambda X: _Aff_eq(X[0], X[1])
    old_eq = lambda X: X[0] == X[1]
    my_int = lambda X: _Aff_intersection(X[0], X[1])
    old_int = lambda X: X[0].intersection(X[1])
    pairs = zip(Aff1, Aff2)

    def TEST(F, pairs):
        start = time.time()
        Results = map(my_eq, pairs)
        end = time.time()
        timing = end - start
        return Results, timing
    
    bv1, t1 = TEST(my_eq, pairs)
    bv2, t2 = TEST(old_eq, pairs)
    if bv1 != bv2: 
        print("Test 1 failed.")
    else:
        print("My eq: %s\nOld eq: %s" % (t1, t2))

    iv1, t3 = TEST(my_int, pairs)
    iv2, t4 = TEST(old_int, pairs)
    if iv1 != iv2: 
        print("Test 2 failed.")
    else:
        print("My int: %s\nOld int: %s" % (t3, t4))
    
    return None


def poset_profile(A):
    import time
    start = time.time()
    P1 = LI.IntersectionPoset(A)
    end = time.time()
    timing1 = end - start
    print("My poset: %s" % (timing1))
    start = time.time()
    P2 = A.intersection_poset()
    end = time.time()
    timing2 = end - start
    print("Old poset: %s" % (timing2))
    assert P1.is_isomorphic(P2)
    if timing1 < timing2:
        print("\t%s percent savings." % (timing2/timing1*100 - 100))
    return P1, P2
    