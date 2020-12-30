# Constructors

We provide a number of constructions of hyperplane arrangements. This uses the default hyperplane arrangement class in SageMath. There is some overlap with the default [library](https://doc.sagemath.org/html/en/reference/discrete_geometry/sage/geometry/hyperplane_arrangement/library.html), and one is encouraged to search there in case your favorite hyperplane arrangement is not in our list. 

Most of our constructions are Coxeter-theoretic&mdash;meaning, they take as input $\mathsf{X}_n$, where $\mathsf{X}$ is the Coxeter type and $n$ is the (Coxeter) rank. 

## CatalanArrangement



## CoxeterArrangement

Input:

- a string or an iterable container of strings.


Output:

- the Coxeter arrangement associated with the strings. 

If just one string is provide, it should be formatted like `'Xn'`, where `X` is a roman letter from $\mathsf{A}$ to $\mathsf{H}$ and `n` is a positive integer. Strings can be separated by one white space like `'Xm Yn'`, and iterable containers of strings need to have strings formatted in this way. 

#### Example (Braid arrangement)

The braid arrangement with $n+1$ hyperplanes is the type $\mathsf{A}_n$ Coxeter arrangement:

\[
    \\{X_i - X_j ~|~ 1\leq i < j \leq n+1\\}.
\]

We construct the braid arrangement with $3$ hyperplanes.

```python
sage: li.CoxeterArrangement("A2")
Arrangement <x1 - x2 | x0 - x1 | x0 - x2>
```

#### Example (Boolean arrangement)

The Boolean arrangement with $n$ hyperplanes is a Coxeter arrangement, isomorphic to $\mathsf{A}_1^n$, the direct sum of $n$ copies of $\mathsf{A}_1$. We construct the rank $5$ Boolean arrangement in two different ways. Note that, when created this way, the ambient dimension is double the rank.

```python
sage: bool5 = "A1 A1 A1 A1 A1"
sage: li.CoxeterArrangement(bool5)
Arrangement of 5 hyperplanes of dimension 10 and rank 5
```

Another construction is as follows.

```python
sage: L = ["A1" for i in range(5)]
sage: L
['A1', 'A1', 'A1', 'A1', 'A1']
sage: li.CoxeterArrangement(L)
Arrangement of 5 hyperplanes of dimension 10 and rank 5
```

#### Example (Coxeter type ${\footnotesize \mathsf{I}_2(m)}$)

In all other Coxeter types, the integer corresponds to the rank. The exception is with type $\mathsf{I}$. If the input is type $\mathsf{I}$, then the integer corresponds to the number of hyperplanes. Here, we give a $\mathbb{Q}$-representation of $\mathsf{I}_2(8)$ as follows.

```python
sage: A = li.CoxeterArrangement("I8")
sage: A
Arrangement of 8 hyperplanes of dimension 2 and rank 2
sage: A.hyperplanes()
(Hyperplane -x0 + 2*x1 + 0,
 Hyperplane -x0 + 3*x1 + 0,
 Hyperplane 0*x0 + x1 + 0,
 Hyperplane x0 - x1 + 0,
 Hyperplane x0 + 0*x1 + 0,
 Hyperplane x0 + x1 + 0,
 Hyperplane x0 + 2*x1 + 0,
 Hyperplane x0 + 3*x1 + 0)
```

## DirectSum

## LinialArrangement 

## PolynomialToArrangement

## ShiArrangement 