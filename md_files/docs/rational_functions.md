# Rational Functions

The main purpose of this package is to explicitly compute the flag Hilbert&ndash;Poincar&#233; series and its specializations like Igusa'a local zeta function. We defer to Maglione&ndash;Voll for the details on the rational functions. 

## BraidArrangementIgusa

**Input**:

- a positive integer.

**Output**:

- the Igusa zeta function associated with the braid arrangement. 

This is a specialized algorithm for the braid arrangement and is significantly faster than `IgusaZetaFunction` on the braid arrangement. This is based off of Lemma 5.14 of Maglione&ndash;Voll.

#### Example (Time comparison)

We compute the Igusa zeta function associated with $\\mathsf{A}_6$ and record the time (on the same machine).
```python
sage: %timeit _ = hi.BraidArrangementIgusa(6)
5.55 ms ± 7.58 μs per loop (mean ± std. dev. of 7 runs, 100 loops each)
sage: %timeit _ = hi.IgusaZetaFunction(hi.CoxeterArrangement("A6"))
4.76 s ± 28 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
```

Thus for $\\mathsf{A}_6$ `BraidArrangementIgusa` is about 1000 times faster than `IgusaZetaFunction`&mdash;of course the latter is also general purpose.

#### Example (Large example)

We compute a factorization of the denominator of the Igusa zeta function associated with $\\mathsf{A}_{10}$.
```python
sage: hi.BraidArrangementIgusa(10).denominator().factor()
(y*t - 1)^5 * (y*t^2 - 1)^2 * (y*t^3 - 1) * (y*t^4 - 1) * (y^2*t^3 - 1)^3 * (y*t^5 - 1) * (y^2*t^4 + y*t^2 + 1)^2 * (y^2*t^5 - 1)^2 * (y^2*t^5 + 1)^2 * (y^2*t^7 - 1) * (y^2*t^9 - 1) * (y^2*t^9 + 1) * (y^2*t^10 + y*t^5 + 1) * (y^2*t^11 - 1) * (y^4*t^12 + y^3*t^9 + y^2*t^6 + y*t^3 + 1) * (y^4*t^14 + y^2*t^7 + 1) * (y^4*t^18 + 1) * (y^6*t^24 + y^5*t^20 + y^4*t^16 + y^3*t^12 + y^2*t^8 + y*t^4 + 1) * (y^6*t^30 + y^3*t^15 + 1) * (y^8*t^44 + y^6*t^33 + y^4*t^22 + y^2*t^11 + 1)
```

Using a program to rewrite the finite geometric progressions in the denominator (like [BRational](https://joshmaglione.com/BRational/)), the denominator takes the following form:
```python
(1 - y*t)^5*(1 - y^2*t^3)^3*(1 - y^3*t^6)^2*(1 - y^4*t^10)^2*(1 - y^5*t^15)*(1 - y^6*t^21)*(1 - y^7*t^28)*(1 - y^8*t^36)*(1 - y^9*t^45)*(1 - y^10*t^55)
```

## CoarseFHPSeries

**Input**:

- `arrangement`: a hyperplane arrangement. Default `None`.
- `matroid`: a matroid. Default `None`.
- `poset`: a poset. Default `None`.
- `R_label`: a function from pairs of elements of a poset into the integers. Default `None`.
- `numerator`: return only the numerator. Default `False`.
- `verbose`: turn on print statements. Default `False`.

**Output**:

- the coarse flag Hilbert&ndash;Poincar&#233; series associated to the graded poset determined by the input. (See the [GradedPoset](graded_posets.md#gradedposet) class.) If `R_label` is given, it is not verified to be an $R$-label, and the algorithm used will be .......

The coarse flag Hilbert&ndash;Poincar&#233; series associated with a graded poset $P$ is defined to be:

\[
    \\mathsf{cfHP}\_{P} (Y, T) 
    = \\sum\_{F\\in\\Delta(P \\setminus \\{\\hat{0}\\})} \\pi\_F(Y) \\left(\\dfrac{T}{1 - T}\\right)^{|F|} 
    = \\dfrac{\\mathcal{N}\_{P}(Y, T)}{(1 - T)^{\\mathrm{rk}(P)}}.
\]

#### Example (Boolean arrangement)

We verify that the Boolean arrangement, $\mathcal{A}$, of rank $n$ satisfies the equation

\[
    \\mathsf{cfHP}\_{\\mathcal{A}}(Y, T) = \\dfrac{(1+ Y)^n E\_n(T)}{(1 - T)^n},
\]

where $E_n(T)$ is the $n$th Eulerian polynomial. We set $n=6$ for this example.

```python 
sage: A = hi.CoxeterArrangement(["A1"]*6)
sage: A
Arrangement of 6 hyperplanes of dimension 12 and rank 6
sage: cfHP = hi.CoarseFHPSeries(A)
sage: cfHP.factor()
(T - 1)^-6 * (T + 1) * (Y + 1)^6 * (T^4 + 56*T^3 + 246*T^2 + 56*T + 1)
```

#### Example (Coxeter at ${\footnotesize Y=1}$)

We verify Theorem D of Maglione--Voll for the Coxeter arrangement of type $\mathsf{D}_5$. Thus, we will show that 

\[
    \\mathsf{cfHP}\_{\\mathsf{D}\_5} (1, T) = 1920\\cdot \\dfrac{1 + 26T + 66T^2 + 26T^3 + T^4}{(1 - T)^5}, 
\]

```python 
sage: A = hi.CoxeterArrangement("D5")
sage: A
Arrangement of 20 hyperplanes of dimension 5 and rank 5
sage: cfHP = hi.CoarseFHPSeries(A)
sage: cfHP.factor()
(-1) * (T - 1)^-5 * (Y + 1) * (Y^4*T^4 + 397*Y^4*T^3 + 19*Y^3*T^4 + 3143*Y^4*T^2 + 3074*Y^3*T^3 + 131*Y^2*T^4 + 3239*Y^4*T + 15624*Y^3*T^2 + 8556*Y^2*T^3 + 389*Y*T^4 + 420*Y^4 + 9694*Y^3*T + 25826*Y^2*T^2 + 9694*Y*T^3 + 420*T^4 + 389*Y^3 + 8556*Y^2*T + 15624*Y*T^2 + 3239*T^3 + 131*Y^2 + 3074*Y*T + 3143*T^2 + 19*Y + 397*T + 1)
```

So we get exactly what we expect:

```python
sage: (cfHP(Y=1)).factor()
(-1920) * (T - 1)^-5 * (T^4 + 26*T^3 + 66*T^2 + 26*T + 1)
```

## FlagHilbertPoincareSeries

**Input**:

- `arrangement`: a hyperplane arrangement. Default `None`.
- `matroid`: a matroid. Default `None`.
- `poset`: a poset. Default `None`.
- `verbose`: turn on print statements. Default `False`.

**Output**:

- the flag Hilbert&ndash;Poincar&#233; series associated to the graded poset determined by the input data. 

The flag Hilbert&ndash;Poincar&#233; series of a graded poset $P$ is defined to be:

\[
    \\mathsf{fHP}\_{P} (Y, \\bm{T}) 
    = \\sum\_{F\\in\\Delta(P \\setminus \\{\\hat{0}\\})} \\pi\_F(Y) \\prod\_{x\\in F} \\frac{T\_x}{1 - T\_x}.
\]

#### Example (Lines through the origin again)

Because of the massive amount of variables in this function, we keep the number of hyperplanes small in this example. We compute the flag Hilbert&ndash;Poincar&#233; series of the arrangement given by 5 lines passes through the origin in a $2$-dimensional vector space. 

```python
sage: K = CyclotomicField(5)
sage: R.<X, Y> = PolynomialRing(K)
sage: f = X**5 - Y**5
sage: A = hi.PolynomialToArrangement(f)
sage: A
Arrangement of 5 hyperplanes of dimension 2 and rank 2
sage: hi.FlagHilbertPoincareSeries(A).factor()
(T6 - 1)^-1 * (T5 - 1)^-1 * (T4 - 1)^-1 * (T3 - 1)^-1 * (T2 - 1)^-1 * (T1 - 1)^-1 * (Y + 1) * (Y*T1*T2*T3*T4*T5 + 4*T1*T2*T3*T4*T5 - Y*T1*T2*T3 - Y*T1*T2*T4 - Y*T1*T3*T4 - Y*T2*T3*T4 - 3*T1*T2*T3*T4 - Y*T1*T2*T5 - Y*T1*T3*T5 - Y*T2*T3*T5 - 3*T1*T2*T3*T5 - Y*T1*T4*T5 - Y*T2*T4*T5 - 3*T1*T2*T4*T5 - Y*T3*T4*T5 - 3*T1*T3*T4*T5 - 3*T2*T3*T4*T5 + 2*Y*T1*T2 + 2*Y*T1*T3 + 2*Y*T2*T3 + 2*T1*T2*T3 + 2*Y*T1*T4 + 2*Y*T2*T4 + 2*T1*T2*T4 + 2*Y*T3*T4 + 2*T1*T3*T4 + 2*T2*T3*T4 + 2*Y*T1*T5 + 2*Y*T2*T5 + 2*T1*T2*T5 + 2*Y*T3*T5 + 2*T1*T3*T5 + 2*T2*T3*T5 + 2*Y*T4*T5 + 2*T1*T4*T5 + 2*T2*T4*T5 + 2*T3*T4*T5 - 3*Y*T1 - 3*Y*T2 - T1*T2 - 3*Y*T3 - T1*T3 - T2*T3 - 3*Y*T4 - T1*T4 - T2*T4 - T3*T4 - 3*Y*T5 - T1*T5 - T2*T5 - T3*T5 - T4*T5 + 4*Y + 1)
```

This is, indeed, equal to 

\[
    \dfrac{1 + Y}{1 - T_6}\left(1 + 4Y + (1 + Y)\sum_{i=1}\dfrac{T_i}{1 - T_i}\right).
\]

## IgusaZetaFunction

**Input**:

- `arrangement`: a hyperplane arrangement. Default `None`.
- `matroid`: a matroid. Default `None`.
- `poset`: a poset. Default `None`.
- `verbose`: turn on print statements. Default `False`.

**Output**:

- the Igusa zeta function associated to the graded poset determined by the input data. 

For a compact discrete valuation ring $\\mathfrak{o}$ and a polynomial $f\\in \\mathfrak{o}[X_1,\dots, X_d]$, Igusa's local zeta function associated with $f$ is 

\[
    Z\_f(s) = \\int\_{\\mathfrak{o}^d} |f(\\bm{X})|^s\\, |\\mathrm{d}\\bm{X}|.
\]

If $Q_\\mathcal{A}$ is the defining polynomial of a hyperplane arrangement, then Igusa's local zeta function associated $\\mathcal{A}$ is $Z_{Q_\\mathcal{A}}(s)$. The output to `IgusaZetaFunction` is a bivariate function in $y = q^{-1}$ and $t = q^{-s}$, where $q$ is the cardinality of the residue field of $\\mathfrak{o}$. This assumes *good reduction*; see Section 1.1 of Maglione&ndash;Voll. 

For data not represented by such a hyperplane arrangement, we define the Igusa zeta function associated with a graded poset $P$ to be 
\[
    Z_P(s) = \\mathsf{fHP}\_P\\left(-q^{-1}, \\left(q^{-g\_x(s)}\\right)\_{x\\in P\\setminus \\{\\hat{0}\\}}\\right) ,
\]
where $g_x(s) = \mathrm{rank}(x) + \\#\\{a\in P \mid a\leqslant x \text{ and $a$ is an atom} \\} \cdot s$. When $P$ is the intersection poset of a hyperplane arrangement $\\mathcal{A}$, then $Z\_P(s)=Z\_{Q\_{\\mathcal{A}}}(s)$, which follows from Theorem B of Maglione&ndash;Voll.

#### Example (Polynomial vs. hyperplane arrangement input)

We demonstrate the two different inputs while showing that polynomials need not have distinct linear factors as is the case with hyperplane arrangements. Let $f(x,y,z) = xy^2z^3$, so that the associated hyperplane arrangement is the Boolean arrangement of rank $3$.

```python 
sage: f = 'x*y^2*z^3'
sage: A = hi.PolynomialToArrangement(f)
sage: A
Arrangement <z | y | x>
```

Now we compare their Igusa zeta functions. The Igusa zeta function associated with $f$ is

```python
sage: Z_f = hi.IgusaZetaFunction(f)
sage: Z_f.factor()
(q - 1)^3/((t^3 - q)*(t^2 - q)*(q - t))
```

which is 

\[
    \dfrac{(1 - q^{-1})^3}{(1 - q^{-1}t) (1 - q^{-1}t^2) (1 - q^{-1}t^3)}.
\]

The Igusa zeta function associated with $\mathcal{A}$ is 

```python
sage: Z_A = hi.IgusaZetaFunction(A)
sage: Z_A.factor()
(q - 1)^3/(q - t)^3
```

which is equal to

\[
    \dfrac{(1 - q^{-1})^3}{(1 - q^{-1}t)^3} .
\]

## TopologicalZetaFuncion

**Input**:

- `arrangement`: a hyperplane arrangement. Default `None`.
- `matroid`: a matroid. Default `None`.
- `poset`: a poset. Default `None`.
- `verbose`: turn on print statements. Default `False`.

**Output**:

- the topological zeta function associated with the graded poset determined by the input data.

For a graded poset $P$, the topological zeta function associated with $P$ is 
\[
    Z\_{P}^{\\mathrm{top}}(s) = \\sum_{F\\in \\Delta(P\\setminus\\{\\hat{0}\\})} \\pi\_{P,F}^\\circ(-1) \\prod_{x\\in F} \\dfrac{1}{g\_x(s)} ,
\]

where $g\_x(s) = \mathrm{rank}(x) + \\#\\{a\in P \mid a\leqslant x \text{ and $a$ is an atom} \\} \cdot s$ and
\[
    \\pi^{\\circ}\_{P, F}(Y) = \\dfrac{\\pi\_F(Y)}{(1 + Y)^{\\# F}}.
\]
When $P$ is the intersection poset associated with a hyperplane arrangement $\\mathcal{A}$, then $Z\_{P}^{\\mathrm{top}}(s) = Z\_{Q_{\\mathcal{A}}}^{\\mathrm{top}}(s)$, which follows from Corollaery 1.5 of Maglione&ndash;Voll.

#### Example (Shi arrangement)

We consider the Shi $\mathsf{A}_2$ arrangement and compute its topological zeta function. The Shi $\mathsf{A}_2$ arrangement is defined to be 
\[
    \mathcal{S} \mathsf{A}_2 = \left\\{X_i - X_j - k ~\middle|~ 1\leqslant i < j\leqslant 3,\; k\in \\{0,1\\}\right\\}.
\]

```python
sage: A = hi.ShiArrangement("A2")
sage: A
Arrangement of 6 hyperplanes of dimension 3 and rank 2
```

The topological zeta function is 

```python
sage: Z = hi.TopologicalZetaFunction(A)
sage: Z.factor()
(s + 1)^-2 * (3*s + 2)^-1 * (12*s^3 + 2*s^2 - 5*s + 2)
```

which is 
\[
    Z\_{\\mathcal{S}\\mathsf{A}\_2}^{\\mathrm{top}}(s) = \\dfrac{2 - 5s + 2s^2 + 12s^3}{(s+1)^2(3s+2)} .
\]
