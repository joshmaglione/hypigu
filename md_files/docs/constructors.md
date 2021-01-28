# Constructors

We provide a number of constructions of hyperplane arrangements. This uses the default hyperplane arrangement class in SageMath. There is some overlap with the default [hyperplane arrangement library](https://doc.sagemath.org/html/en/reference/discrete_geometry/sage/geometry/hyperplane_arrangement/library.html), and one is encouraged to search there in case your favorite hyperplane arrangement is not in our list. 

Most of our constructions are Coxeter-theoretic&mdash;meaning, they take as input $\mathsf{X}$ and $n$, where $\mathsf{X}$ is the Coxeter type and $n$ is the (Coxeter) rank. For example the Catalan arrangement of type $\mathsf{X}_n$ can be defined for all Coxeter arrangements like its defined for the (type $\mathsf{A}$) braid arrangement. 

## CatalanArrangement

**Input**:

- a string or an iterable container of strings.


**Output**:

- the associated Catalan arrangement. 

If just one string is provided, it should be formatted like `'Xn'`, where `X` is a roman letter from $\mathsf{A}$ to $\mathsf{H}$ and `n` is a positive integer. Iterable containers of strings need to have strings formatted in this way. Strings can also be separated by one white space like `'Xm Yn'` instead of being in an iterable container. 

If $\mathcal{A}$ is a Coxeter arrangement of type $\mathsf{X}_n$, then the *Catalan arrangement of type $\mathsf{X}_n$* is 

\[ 
    \mathcal{C} = \\{L - 1 ~|~ L\in\mathcal{A} \\} \cup \mathcal{A} \cup \\{L + 1 ~|~ L\in\mathcal{A} \\}.
\]

#### Example (*The* Catalan arrangement)

The Catalan arrangement is defined as 

\[
    \\{X_i - X_j + \varepsilon ~|~ 1\leq i < j \leq n+1,\; \varepsilon\in\\{-1, 0, 1\\} \\}.
\]

We can quickly construct this as a type $\mathsf{A}_n$ Catalan arrangement for $n=2$. 

```python
sage: A = hi.CatalanArrangement("A2")
sage: A
Arrangement of 9 hyperplanes of dimension 3 and rank 2
sage: A.hyperplanes()
(Hyperplane 0*x0 + x1 - x2 - 1,
 Hyperplane 0*x0 + x1 - x2 + 0,
 Hyperplane 0*x0 + x1 - x2 + 1,
 Hyperplane x0 - x1 + 0*x2 - 1,
 Hyperplane x0 - x1 + 0*x2 + 0,
 Hyperplane x0 - x1 + 0*x2 + 1,
 Hyperplane x0 + 0*x1 - x2 - 1,
 Hyperplane x0 + 0*x1 - x2 + 0,
 Hyperplane x0 + 0*x1 - x2 + 1)
```

#### Example (${\footnotesize \mathsf{B}_4}$-Catalan arrangement)

The Coxeter arrangement of type $\mathsf{B}_4$ has $16$ hyperplanes in $\mathbb{Q}^4$. Thus, the associated Catalan arrangement is non-central with $48$ hyperplanes in $\mathbb{Q}^4$. We verify this.

```python
sage: A = hi.CatalanArrangement("B4")
sage: A
Arrangement of 48 hyperplanes of dimension 4 and rank 4
sage: A.is_central()
False
```

## CoxeterArrangement

**Input**:

- a string or an iterable container of strings.

**Output**:

- the associated Coxeter arrangement. 

If just one string is provided, it should be formatted like `'Xn'`, where `X` is a roman letter from $\mathsf{A}$ to $\mathsf{H}$ and `n` is a positive integer. Iterable containers of strings need to have strings formatted in this way. Strings can also be separated by one white space like `'Xm Yn'` instead of being in an iterable container. 

#### Example (Braid arrangement)

The braid arrangement with $n+1$ hyperplanes is equivalent to the type $\mathsf{A}_n$ Coxeter arrangement:

\[
    \\{X_i - X_j ~|~ 1\leq i < j \leq n+1\\}.
\]

We construct the braid arrangement with $3$ hyperplanes.

```python
sage: hi.CoxeterArrangement("A2")
Arrangement <x1 - x2 | x0 - x1 | x0 - x2>
```

#### Example (Boolean arrangement)

One should consider the other two constructions of the Boolean arrangement via [direct sums](#example-boolean-arrangement-again) and [polynomials](#example-boolean-arrangement-yet-again). The Boolean arrangement with $n$ hyperplanes is a Coxeter arrangement, equivalent to $\mathsf{A}_1^n$, the direct sum of $n$ copies of $\mathsf{A}_1$. We construct the rank $5$ Boolean arrangement in two different ways. Note that, when created this way, the ambient dimension is double the rank.

```python
sage: bool5 = "A1 A1 A1 A1 A1"
sage: hi.CoxeterArrangement(bool5)
Arrangement of 5 hyperplanes of dimension 10 and rank 5
```

Another construction is as follows.

```python
sage: L = ["A1" for i in range(5)]
sage: L
['A1', 'A1', 'A1', 'A1', 'A1']
sage: hi.CoxeterArrangement(L)
Arrangement of 5 hyperplanes of dimension 10 and rank 5
```

#### Example (Coxeter type ${\footnotesize \mathsf{I}_2(m)}$)

In all other Coxeter types, the integer corresponds to the rank. The exception is with type $\mathsf{I}$. If the input is type $\mathsf{I}$, then the integer corresponds to the number of hyperplanes. Here, we give a $\mathbb{Q}$-representation of $\mathsf{I}_2(8)$ as follows.

```python
sage: A = hi.CoxeterArrangement("I8")
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

**Input**:

- an iterable container of hyperplane arrangements.

**Output**:

- the direct sum arrangement. 

The direct sum arrangement is also known as the product arrangement. 

#### Example (Boolean arrangement again)

One should compare this example with the other Boolean examples: as a [Coxeter arrangement](#example-boolean-arrangement) and as an arrangement from a [polynomial](#example-boolean-arrangement-yet-again). First we create the irreducible factor: the origin on the rational line. 

```python 
sage: H = HyperplaneArrangements(QQ, 'x')
sage: A = H([0, 1])
sage: A
Arrangement <x>
```

Now we create higher rank Boolean arrangements in two slightly different ways. The rank $4$ Boolean arrangement is given as follows.

```python
sage: hi.DirectSum(A, A, A, A)
Arrangement <x3 | x2 | x1 | x0>
```

And the rank $256$ Boolean arrangement is 

```python 
sage: L = [A for i in range(256)]
sage: hi.DirectSum(L)
Arrangement of 256 hyperplanes of dimension 256 and rank 256
```

## LinialArrangement 

**Input**:

- a string or an iterable container of strings.

**Output**:

- the associated Linial arrangement. 

If just one string is provided, it should be formatted like `'Xn'`, where `X` is a roman letter from $\mathsf{A}$ to $\mathsf{H}$ and `n` is a positive integer. Iterable containers of strings need to have strings formatted in this way. Strings can also be separated by one white space like `'Xm Yn'` instead of being in an iterable container. 

If $\mathcal{A}$ is a Coxeter arrangement of type $\mathsf{X}_n$, then the *Linial arrangement of type $\mathsf{X}_n$* is 

\[ 
    \mathcal{L} = \\{L - 1 ~|~ L\in\mathcal{A} \\} .
\]

#### Example (*The* Linial arrangement)

Usually the Linial arrangement is defined without reference to a Coxeter type and is given as the type-$\mathsf{A}$ version above. We construct the Linial arrangement of type $\mathsf{A}_2$.

```pyhton
sage: A = hi.LinialArrangement("A2")
sage: A
Arrangement <x1 - x2 + 1 | x0 - x1 + 1 | x0 - x2 + 1>
```

## PolynomialToArrangement

**Input**:

- a polynomial or symbolic expression.

**Output**:

- the hyperplane arrangement associated with the linear factors of the given polynomial.

We require that the given polynomial *only* have linear factors. Strings are also acceptable input and will be interpreted by SageMath (so they should be formatted accordingly). The underlying field for such symbolic expressions is assumed to be $\mathbb{Q}$. 

#### Example (Boolean arrangement yet again)

Compare this construction with the Boolean construction as a [Coxeter arrangement](#example-boolean-arrangement) and as a [direct sum](#example-boolean-arrangement-again). The Boolean arrangement is equivalent to the arrangement of coordinate hyperplanes, so the corresponding polynomial is a monomial. We construct the rank $5$ Boolean arrangement.

```python
sage: a, b, c, d, e = var('a b c d e')
sage: f = a*b*c*d*e
sage: A = hi.PolynomialToArrangement(f)
sage: A
Arrangement of 5 hyperplanes of dimension 5 and rank 5
```

Equivalently, one can input just the string.

```python
sage: B = hi.PolynomialToArrangement('a*b*c*d*e')
sage: A == B
True
```

#### Example (Cyclotomic)

Let $\zeta_5$ be a primitive $5$th root of unity. We construct the arrangement of $5$ lines through the origin with underlying field $\mathbb{Q}(\zeta_5)$. This can be accomplished by constructing the hyperplane arrangement associated to the polynomial $X^5 - Y^5\in \mathbb{Q}(\zeta_5)[X, Y]$. Note that this arrangement is equivalent to the Coxeter arrangement of type $\mathsf{I}_2(5)\cong \mathsf{H}_2$; see the [$\mathsf{I}_2(m)$ Example](#example-coxeter-type-footnotesize-mathsfi_2m). 

```python 
sage: K = CyclotomicField(5)
sage: R.<X, Y> = PolynomialRing(K)
sage: f = X**5 - Y**5
sage: A = hi.PolynomialToArrangement(f)
sage: A
Arrangement of 5 hyperplanes of dimension 2 and rank 2
sage: A.hyperplanes()
(Hyperplane X + (-zeta5^2)*Y + 0,
 Hyperplane X + (zeta5^3 + zeta5^2 + zeta5 + 1)*Y + 0,
 Hyperplane X + (-zeta5^3)*Y + 0,
 Hyperplane X + (-zeta5)*Y + 0,
 Hyperplane X + (-1)*Y + 0)
```

Note that if the underlying field were just, say, $\mathbb{Q}$, then this would result in an error since $f$ is not a product of linear factors. 

## ShiArrangement 

**Input**:

- a string or an iterable container of strings.

**Output**:

- the associated Shi arrangement. 

If just one string is provided, it should be formatted like `'Xn'`, where `X` is a roman letter from $\mathsf{A}$ to $\mathsf{H}$ and `n` is a positive integer. Iterable containers of strings need to have strings formatted in this way. Strings can also be separated by one white space like `'Xm Yn'` instead of being in an iterable container. 

If $\mathcal{A}$ is a Coxeter arrangement of type $\mathsf{X}_n$, then the *Shi arrangement of type $\mathsf{X}_n$* is 

\[ 
    \mathcal{S} = \mathcal{A} \cup \\{L - 1 ~|~ L\in\mathcal{A} \\} .
\]

#### Example (*The* Shi arrangement)

Like with some of our other constructors, the usual definition makes no reference to Coxeter types, so the Shi arrangement is equal to the type-$\mathsf{A}$ Shi arrangement defined above. We can easily construct this for $\mathsf{A}_2$.

```python
sage: A = hi.ShiArrangement("A2")
sage: A
Arrangement of 6 hyperplanes of dimension 3 and rank 2
```