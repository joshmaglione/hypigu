# Graded Posets

We include a class called `GradedPoset`. This class keeps track of basic information we need to compute the various generating functions associated with a hyperplane arrangement.

## GradedPoset

**Input**:

- `arrangement`: a hyperplane arrangement. Default `None`.
- `matroid`: a matroid. Default `None`.
- `poset`: a graded poset with a bottom element. Default `None`.
- `R_label`: a function from the covering relations (pairs of `poset`) to the integers. Default `None`.

**Output**: 

- the corresponding graded poset. If an arrangement is given, the poset comes from the intersection poset. If a matroid is given, the poset comes from the lattice of flats. $R$-labels are not checked, so we assume whatever is given is correct.

### Attributes 

The class `GradedPoset` has four attributes, which are the same as the keyword input.

#### Example (Lattice of braid arrangement)

We construct the lattice of flats for the braid arrangement in $\mathbb{R}^4$.

```python
sage: A = hi.CoxeterArrangement("A3")
sage: A
Arrangement of 6 hyperplanes of dimension 4 and rank 3
sage: GP = hi.GradedPoset(arrangement=A)
sage: GP
An R-labeled graded poset with 15 elements built from the intersection poset of
Arrangement of 6 hyperplanes of dimension 4 and rank 3
```

Now we look at the data stored in the attributes. We display the hyperplanes in the arrangement.

```python
sage: GP.arrangement.hyperplanes()
(Hyperplane 0*x0 + 0*x1 + x2 - x3 + 0,
 Hyperplane 0*x0 + x1 - x2 + 0*x3 + 0,
 Hyperplane 0*x0 + x1 + 0*x2 - x3 + 0,
 Hyperplane x0 - x1 + 0*x2 + 0*x3 + 0,
 Hyperplane x0 + 0*x1 - x2 + 0*x3 + 0,
 Hyperplane x0 + 0*x1 + 0*x2 - x3 + 0)
```

We display the poset as an image.

```python
sage: GP.poset
Finite lattice containing 15 elements
```

![](A3.png)


We see that hyperplanes 0, 4, and 5 intersect in a codimension $2$ subspace. 

```python
sage: [GP.arrangement.hyperplanes()[i] for i in [0,4,5]]
[Hyperplane 0*x0 + 0*x1 + x2 - x3 + 0,
 Hyperplane x0 + 0*x1 - x2 + 0*x3 + 0,
 Hyperplane x0 + 0*x1 + 0*x2 - x3 + 0]
```

## .atoms

**Output**:

- the atoms of the underlying poset. 

## .interval

**Input**: 

- `bottom`: the bottom element of the desired interval. Default `None`.
- `top`: the top element of the desired interval. Default `None`.

**Output**:

- the closed interval from `bottom` (if `None` the unique bottom element is used) to `top` (if `None` all suitable maximal elements are used).

## .Poincare_polynomial

**Output**:

- the Poincar&#233; polynomial of the graded poset. This is defined to be

\[
	\\pi\_P(Y) = \\sum\_{x\\in P} |\\mu(\\hat{0}, x)|\\cdot Y^{\\mathrm{rank}(x)}.
\]

## .show

No output given. This displays the underlying intersection poset.
