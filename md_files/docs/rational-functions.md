# Rational Functions

The main purpose of this package is to explicitly compute the flag Hilbert&ndash;Poincar&#233; series and its specializations like Igusa'a local zeta function. We keep variable letters consistent with Maglione&ndash;Voll; the exception is that we replace $q^{-s_x}$ by $t_x$, for some label $x$. We define all of the rational functions, but defer to Maglione&ndash;Voll for the details. 

We omit the analogous definitions for matroids. Since all of the rational functions here are related to the flag Hilbert&ndash;Poincar&#233; series, the matroid versions use the lattice of flats of the given matroid and possibly Theorem B of Maglione&ndash;Voll.

All of these functions contains many *optional* parameters. With the exception of [TopologicalZetaFunction](#topologicalzetafuncion), these are present to either provide print statements or save on computation time by using data previously computed. Unless these data have been computed, one should leave such parameters set to `None`.

For example, if one wants to compute the combinatorial skeleton, the Igusa zeta function, and the topological zeta function of $\mathcal{A}$, then it would be beneficial to first compute the [lattice of flats](https://joshmaglione.github.io/hypigu/lattices/#latticeofflats) and pass it to the three rational functions. This way, the lattice of flats is computed only once.

```python
sage: A = hi.CoxeterArrangement("A5")
sage: L = hi.LatticeOfFlats(A)
sage: CS = hi.CombinatorialSkeleton(A, lattice_of_flats=L)
sage: Z = hi.IgusaZetaFunction(A, lattice_of_flats=L)
sage: TZ = hi.TopologicalZetaFunction(A, lattice_of_flats=L)
```

## AnalyticZetaFunction

**Input**:

- a hyperplane arrangement $\mathcal{A}$,
- `matroid=None` : a matroid, 
- `lattice_of_flats=None` : the lattice of flats of $\mathcal{A}$,
- `int_poset=None` : the intersection poset of $\mathcal{A}$,
- `verbose=False` : turn on print statements.

**Output**:

- the analytic zeta function associated to $\mathcal{A}$. 

Given a suitable $\mathfrak{o}$-representation of a $d$-dimensional hyperplane arrangement $\mathcal{A}$, where $\mathfrak{o}$ is a compact discrete valuation ring, the analytic zeta function is defined to be the integral:

\[
    \zeta_{\mathcal{A}(\mathfrak{o})}(\bm{s}) = \int_{\mathfrak{o}^d} \prod_{x\in\widetilde{\mathcal{L}}(\mathcal{A})} \\| \mathcal{A}_x\\|^{s_x} \\, |\mathrm{d}\bm{X}|. 
\]

The parameter `lattice_of_flats` can be used to give the lattice of flats of $\mathcal{A}$, computed by [LatticeOfFlats](https://joshmaglione.github.io/hypigu/lattices/#latticeofflats); otherwise this parameter should stay set to `None`. The paramater `int_poset` can be used to give in the intersection poset of $\mathcal{A}$; otherwise this parameter should stay set to `None`.

#### Example (Lines through the origin)

We compute one of the analytic zeta functions in Section 4.2 of Maglione&ndash;Voll. Let $m\geq 2$ and $\zeta_m$ a primitive $m$th root of unity. We define an arrangement of $m$ lines through the origin, given as the linear factors of $X^m-Y^m$ in $\mathbb{Q}(\zeta_m)$. We set $m=5$ for this example. 

```python
sage: K = CyclotomicField(5)
sage: R.<X, Y> = PolynomialRing(K)
sage: f = X**5 - Y**5
sage: A = hi.PolynomialToArrangement(f)
sage: A
Arrangement of 5 hyperplanes of dimension 2 and rank 2
```

The analytic zeta function is then

```python 
sage: Z = hi.AnalyticZetaFunction(A)
sage: Z
-((4/q - 1)*(1/q - 1) - t1*(1/q - 1)^2/(q*(t1/q - 1)) - t2*(1/q - 1)^2/(q*(t2/q - 1)) - t3*(1/q - 1)^2/(q*(t3/q - 1)) - t4*(1/q - 1)^2/(q*(t4/q - 1)) - t5*(1/q - 1)^2/(q*(t5/q - 1)))/(t1*t2*t3*t4*t5*t6/q^2 - 1)
```

which is indeed

\[
    \dfrac{1 - q^{-1}}{1 - q^{-2-s_{\hat{1}} - s_1-\cdots -s_5}} \left(1 - 4q^{-1} + (1 - q^{-1}) \sum_{i=1}^5 \dfrac{q^{-1-s_i}}{1 - q^{-1-s_i}}\right) .
\]

## AtomZetaFunction

**Input**:

- a hyperplane arrangement $\mathcal{A}$,
- `matroid=None` : a matroid, 
- `lattice_of_flats=None` : the lattice of flats of $\mathcal{A}$,
- `int_poset=None` : the intersection poset of $\mathcal{A}$,
- `verbose=False` : turn on print statements.

**Output**:

- the atom zeta function associated to $\mathcal{A}$. 

Given a suitable $\mathfrak{o}$-representation of a $d$-dimensional hyperplane arrangement $\mathcal{A}$, where $\mathfrak{o}$ is a compact discrete valuation ring, the atom zeta function is defined to be the integral:

\[
    \zeta_{\mathcal{A}(\mathfrak{o})}^{\mathrm{at}}(\textbf{s}) = \int_{\mathfrak{o}^d} \prod_{L\in\mathcal{A}(\mathfrak{o})} |L(\bm{X})|^{s_L} \\, |\mathrm{d}\bm{X}|. 
\]

The parameter `lattice_of_flats` can be used to give the lattice of flats of $\mathcal{A}$, computed by [LatticeOfFlats](https://joshmaglione.github.io/hypigu/lattices/#latticeofflats); otherwise this parameter should stay set to `None`. The paramater `int_poset` can be used to give in the intersection poset of $\mathcal{A}$; otherwise this parameter should stay set to `None`.

#### Example (Atom zeta function for braid arrangement)

We compute the atom zeta function for the braid arrangement, $\mathcal{A}$, in $\mathbb{R}^4$, which has $6$ hyperplanes, so $\zeta_{\mathcal{A}(\mathfrak{o})}^{\mathrm{at}}(\bm{s})$ has $6$ variables. First we construct the braid arrangement as a Coxeter arrangement of type $\mathsf{A}_3$. 

```python
sage: A = hi.CoxeterArrangement("A3")
sage: A
Arrangement of 6 hyperplanes of dimension 4 and rank 3
```

Now we construct the atom zeta function of $\mathcal{A}$.

```python
sage: Z = hi.AtomZetaFunction(A)
sage: Z
-(((2/q - 1)*(1/q - 1) - t1*(1/q - 1)^2/(q*(t1/q - 1)) - t2*(1/q - 1)^2/(q*(t2/q - 1)) - t3*(1/q - 1)^2/(q*(t3/q - 1)))*t1*t2*t3*(1/q - 1)/(q^2*(t1*t2*t3/q^2 - 1)) + ((2/q - 1)*(1/q - 1) - t2*(1/q - 1)^2/(q*(t2/q - 1)) - t4*(1/q - 1)^2/(q*(t4/q - 1)) - t5*(1/q - 1)^2/(q*(t5/q - 1)))*t2*t4*t5*(1/q - 1)/(q^2*(t2*t4*t5/q^2 - 1)) + ((2/q - 1)*(1/q - 1) - t3*(1/q - 1)^2/(q*(t3/q - 1)) - t4*(1/q - 1)^2/(q*(t4/q - 1)) - t6*(1/q - 1)^2/(q*(t6/q - 1)))*t3*t4*t6*(1/q - 1)/(q^2*(t3*t4*t6/q^2 - 1)) + ((2/q - 1)*(1/q - 1) - t1*(1/q - 1)^2/(q*(t1/q - 1)) - t5*(1/q - 1)^2/(q*(t5/q - 1)) - t6*(1/q - 1)^2/(q*(t6/q - 1)))*t1*t5*t6*(1/q - 1)/(q^2*(t1*t5*t6/q^2 - 1)) + ((1/q - 1)^2 - t1*(1/q - 1)^2/(q*(t1/q - 1)) - t4*(1/q - 1)^2/(q*(t4/q - 1)))*t1*t4*(1/q - 1)/(q^2*(t1*t4/q^2 - 1)) + ((1/q - 1)^2 - t3*(1/q - 1)^2/(q*(t3/q - 1)) - t5*(1/q - 1)^2/(q*(t5/q - 1)))*t3*t5*(1/q - 1)/(q^2*(t3*t5/q^2 - 1)) + ((1/q - 1)^2 - t2*(1/q - 1)^2/(q*(t2/q - 1)) - t6*(1/q - 1)^2/(q*(t6/q - 1)))*t2*t6*(1/q - 1)/(q^2*(t2*t6/q^2 - 1)) - t1*(3/q - 2/q^2 - 1)*(1/q - 1)/(q*(t1/q - 1)) - t2*(3/q - 2/q^2 - 1)*(1/q - 1)/(q*(t2/q - 1)) - t3*(3/q - 2/q^2 - 1)*(1/q - 1)/(q*(t3/q - 1)) - t4*(3/q - 2/q^2 - 1)*(1/q - 1)/(q*(t4/q - 1)) - t5*(3/q - 2/q^2 - 1)*(1/q - 1)/(q*(t5/q - 1)) - t6*(3/q - 2/q^2 - 1)*(1/q - 1)/(q*(t6/q - 1)) - 6/q + 11/q^2 - 6/q^3 + 1)/(t1*t2*t3*t4*t5*t6/q^3 - 1)
```

Expressing $\zeta_{\mathcal{A}(\mathfrak{o})}^{\mathrm{at}}(\bm{s})$ as a quotient of polynomials requires too much text space for this example, so we will just write the denominator:

```python
sage: Z.numerator_denominator()[1]
-(t1*t2*t3*t4*t5*t6 - q^3)*(t1*t2*t3 - q^2)*(t2*t4*t5 - q^2)*(t3*t4*t6 - q^2)*(t1*t5*t6 - q^2)*(q - t1)*(q - t2)*(q - t3)*(q - t4)*(q - t5)*(q - t6)
```

which, up to multiples of $q$, is

\[
    \left(1-q^{-3}t_1\cdots t_6\right) \left(1-q^{-2}t_1t_2t_3\right) \left(1-q^{-2}t_2t_4t_5\right) \left(1-q^{-2}t_3t_4t_6\right) \left(1-q^{-2}t_1t_5t_6\right) \prod_{i=1}^6 \left(1 - q^{-1}t_i\right).
\]


#### Example (Multiplicativity of atom zeta functions)

We verify that the atom zeta function is multiplicative, which follows from Fubini's theorem. We demonstrate this using the Boolean arrangement: $\mathsf{A}_1^n$. 

```python
sage: A = hi.CoxeterArrangement("A1")
sage: A
Arrangement <x0 - x1>
sage: B = hi.DirectSum([A for i in range(5)])
sage: B
Arrangement of 5 hyperplanes of dimension 10 and rank 5
```

Now we compute the atom zeta functions of $\mathcal{A}$ and $\mathcal{B}$.

```python 
sage: Z_A = hi.AtomZetaFunction(A)
sage: Z_B = hi.AtomZetaFunction(B)
```

The first atom zeta function is simple:

```python
sage: Z_A
(1/q - 1)/(t1/q - 1)
```

and equal to 

\[
    \dfrac{1-q^{-1}}{1-q^{-1}t_1}.
\]

The second atom zeta function is 

```python
sage: Z_B.numerator_denominator()
(q^5 - 5*q^4 + 10*q^3 - 10*q^2 + 5*q - 1,
 (q - t1)*(q - t2)*(q - t3)*(q - t4)*(q - t5))
```

which is equal to 

\[
    \prod_{i=1}^5\dfrac{1-q^{-1}}{1 - q^{-1}t_i}.
\]

## CoarseFlagHPSeries

**Input**:

- a hyperplane arrangement $\mathcal{A}$,
- `matroid=None` : a matroid, 
- `lattice_of_flats=None` : the lattice of flats of $\mathcal{A}$,
- `int_poset=None` : the intersection poset of $\mathcal{A}$,
- `numerator=False` : only return the numerator $\mathcal{N}_{\mathcal{A}}(Y, T)$,
- `verbose=False` : turn on print statements.

**Output**:

- the coarse flag Hilbert&ndash;Poincar&#233; series associated to $\mathcal{A}$. 

The coarse flag Hilbert&ndash;Poincar&#233; series of $\mathcal{A}$ is defined to be:

\[
    cfHP_{\mathcal{A}} (Y, T) 
    = \sum_{F\in\Delta(\widetilde{\mathcal{L}}(\mathcal{A}))} \pi_F(Y) \left(\dfrac{T}{1 - T}\right)^{|F|} 
    = \dfrac{\mathcal{N}_{\mathcal{A}}(Y, T)}{(1 - T)^{\mathrm{rk}(\mathcal{A})}}.
\]

The parameter `lattice_of_flats` can be used to give the lattice of flats of $\mathcal{A}$, computed by [LatticeOfFlats](https://joshmaglione.github.io/hypigu/lattices/#latticeofflats); otherwise this parameter should stay set to `None`. The paramater `int_poset` can be used to give in the intersection poset of $\mathcal{A}$; otherwise this parameter should stay set to `None`.

#### Example (Boolean skeleton)

We verify that the Boolean arrangement, $\mathcal{A}$, of rank $n$ satisfies the equation

\[
    {\sf cfHP}_{\mathcal{A}}(Y, T) = \dfrac{(1+ Y)^nE_n(T)}{(1 - T)^n},
\]

where $E_n(T)$ is the $n$th Eulerian polynomial. We set $n=6$ for this example.

```python 
sage: A = hi.CoxeterArrangement(["A1" for i in range(6)])
sage: A
Arrangement of 6 hyperplanes of dimension 12 and rank 6
sage: S = hi.CoarseFlagHPSeries(A)
sage: S.factor()
(T^4 + 56*T^3 + 246*T^2 + 56*T + 1)*(T + 1)*(Y + 1)^6/(T - 1)^6
```

#### Example (Coxeter skeletons at ${\footnotesize Y=1}$)

We verify Theorem D of Maglione--Voll for the Coxeter arrangement of type $\mathsf{D}_5$. Thus, we will show that 

\[
    {\sf cfHP}_{\mathsf{D}_5} (1, T) = 1920\cdot \dfrac{1 + 26T + 66T^2 + 26T^3 + T^4}{(1 - T)^5}, 
\]

```python 
sage: A = hi.CoxeterArrangement("D5")
sage: A
Arrangement of 20 hyperplanes of dimension 5 and rank 5
sage: S = hi.CoarseFlagHPSeries(A)
sage: S.factor()
-(T^4*Y^4 + 19*T^4*Y^3 + 397*T^3*Y^4 + 131*T^4*Y^2 + 3074*T^3*Y^3 + 3143*T^2*Y^4 + 389*T^4*Y + 8556*T^3*Y^2 + 15624*T^2*Y^3 + 3239*T*Y^4 + 420*T^4 + 9694*T^3*Y + 25826*T^2*Y^2 + 9694*T*Y^3 + 420*Y^4 + 3239*T^3 + 15624*T^2*Y + 8556*T*Y^2 + 389*Y^3 + 3143*T^2 + 3074*T*Y + 131*Y^2 + 397*T + 19*Y + 1)*(Y + 1)/(T - 1)^5
```

So we get exactly what we expect:

```python
sage: S(Y=1).factor()/1920
-(T^4 + 26*T^3 + 66*T^2 + 26*T + 1)/(T - 1)^5
```

## FlagHilbertPoincareSeries

**Input**:

- a hyperplane arrangement $\mathcal{A}$,
- `matroid=None` : a matroid, 
- `lattice_of_flats=None` : the lattice of flats of $\mathcal{A}$,
- `int_poset=None` : the intersection poset of $\mathcal{A}$,
- `verbose=False` : turn on print statements.

**Output**:

- the flag Hilbert&ndash;Poincar&#233; series associated to $\mathcal{A}$. 

The flag Hilbert&ndash;Poincar&#233; series of $\mathcal{A}$ is defined to be:

\[
    fHP_{\mathcal{A}} (Y, \bm{T}) 
    = \sum_{F\in\Delta(\widetilde{\mathcal{L}}(\mathcal{A}))} \pi_F(Y) \prod_{x\in F} \frac{T_x}{1 - T_x}.
\]

The parameter `lattice_of_flats` can be used to give the lattice of flats of $\mathcal{A}$, computed by [LatticeOfFlats](https://joshmaglione.github.io/hypigu/lattices/#latticeofflats); otherwise this parameter should stay set to `None`. The paramater `int_poset` can be used to give in the intersection poset of $\mathcal{A}$; otherwise this parameter should stay set to `None`.

#### Example (Lines through the origin again)

Because of the massive amount of variables in this function (and the analytic zeta function), we keep the number of hyperplanes small in this example. We compute the flag Hilbert&ndash;Poincar&#233; series of the same arrangement given in the analytic zeta function [example](#example-lines-through-the-origin), so we will not redo the construction of $\mathcal{A}$. 

```python
sage: A
Arrangement of 5 hyperplanes of dimension 2 and rank 2
sage: hi.FlagHilbertPoincareSeries(A)
-((4*Y + 1)*(Y + 1) - T1*(Y + 1)^2/(T1 - 1) - T2*(Y + 1)^2/(T2 - 1) - T3*(Y + 1)^2/(T3 - 1) - T4*(Y + 1)^2/(T4 - 1) - T5*(Y + 1)^2/(T5 - 1))/(T6 - 1)
```

This is, indeed, equal to 

\[
    \dfrac{1 + Y}{1 - T_6}\left(1 + 4Y + (1 + Y)\sum_{i=1}\dfrac{T_i}{1 - T_i}\right).
\]

## IgusaZetaFunction

**Input**:

- a hyperplane arrangement $\mathcal{A}$ or a polynomial $f$,
- `matroid=None` : a matroid, 
- `lattice_of_flats=None` : the lattice of flats of $\mathcal{A}$,
- `int_poset=None` : the intersection poset of $\mathcal{A}$,
- `verbose=False` : turn on print statements.

**Output**:

- Igusa's local zeta function associated to either $\mathcal{A}$ or $f$. 

If a polynomial, $f$, is given, we require that $f$ be the product of linear factors. Symbolic expressions and strings are fine as well, provided SageMath interprets them as a polynomial. This kind of input should be acceptable for [PolynomialToArrangement](https://joshmaglione.github.io/hypigu/constructors/#polynomialtoarrangement).

For a compact discrete valuation ring $\mathfrak{o}$ and a polynomial $f\in \mathfrak{o}[X_1,\dots, X_d]$, Igusa's local zeta function associated with $f$ is 

\[
    Z_f(s) = \int_{\mathfrak{o}^d} |f(\bm{X})|^s\, |\mathrm{d}\bm{X}|.
\]

If $Q_\mathcal{A}$ is the defining polynomial of a hyperplane arrangement, then Igusa's local zeta function associated $\mathcal{A}$ is $Z_{Q_\mathcal{A}}(s)$. 

The parameter `lattice_of_flats` can be used to give the lattice of flats of $\mathcal{A}$, computed by [LatticeOfFlats](https://joshmaglione.github.io/hypigu/lattices/#latticeofflats); otherwise this parameter should stay set to `None`. The paramater `int_poset` can be used to give in the intersection poset of $\mathcal{A}$; otherwise this parameter should stay set to `None`.

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

- a hyperplane arrangement $\mathcal{A}$ or a polynomial $f$,
- `matroid=None` : a matroid, 
- `multivariate=False` : return the *multivariate* zeta function associated with $\mathcal{A}$,
- `atom=False` : return the *atom specialization* of the multivariate zeta function associated with $\mathcal{A}$,
- `lattice_of_flats=None` : the lattice of flats of $\mathcal{A}$,
- `int_poset=None` : the intersection poset of $\mathcal{A}$,
- `verbose=False` : turn on print statements.

**Output**:

- the topological zeta function associated to either $\mathcal{A}$ or $f$. 

If a polynomial, $f$, is given, we require that $f$ be the product of linear factors. Symbolic expressions and strings are fine as well, provided SageMath interprets them as a polynomial. This kind of input should be acceptable for [PolynomialToArrangement](https://joshmaglione.github.io/hypigu/constructors/#polynomialtoarrangement).

For a hyperplane arrangement $\mathcal{A}$, the multivariate topological zeta function associated with $\mathcal{A}$ is 
\[
    \zeta_{\mathcal{A}}^{\mathrm{top}}(\bm{s}) = \sum_{F\in \Delta(\widetilde{\mathcal{L}}(\mathcal{A}))} \pi_{\mathcal{A},F}^\circ(-1) \prod_{x\in F} \dfrac{1}{\mathrm{rk}(x) + \sum_{y\in\widetilde{\mathcal{L}}(\mathcal{A}_x)}s_y}  .
\]

Depending on the parameters `multivariate` and `atom`, different topological zeta functions are returned. If `multivariate=True` and `atom=False`, then the multivariate topological zeta function is returned. If `mutlivariate=True` and `atom=True`, then the *atom specialization* is returned; namely,
\[
    \zeta_{\mathcal{A}}^{\mathrm{top},\mathrm{at}}(\bm{s}) = 
    \zeta_{\mathcal{A}}^{\mathrm{top}}\left((s_x\cdot \delta_{|A_x|=1})_{x\in\widetilde{\mathcal{L}}(\mathcal{A})}\right).
\]

Lastly, if `multivariate=False`, then the (univariate) topological zeta function is returned, which is defined to be
\[
    Z_{\mathcal{A}}^{\mathrm{top}}(s) = \zeta_{\mathcal{A}}^{\mathrm{top},\mathrm{at}}\left((s)_{L\in\mathcal{A}}\right) .
\]

The parameter `lattice_of_flats` can be used to give the lattice of flats of $\mathcal{A}$, computed by [LatticeOfFlats](https://joshmaglione.github.io/hypigu/lattices/#latticeofflats); otherwise this parameter should stay set to `None`. The paramater `int_poset` can be used to give in the intersection poset of $\mathcal{A}$; otherwise this parameter should stay set to `None`.

#### Example (Shi arrangement)

We consider the Shi $\mathsf{A}_2$ arrangement and compute its topological zeta function. The Shi $\mathsf{A}_2$ arrangement is defined to be 
\[
    \mathcal{S} \mathsf{A}_2 = \left\\{X_i - X_j - k ~\middle|~ 1\leq i < j\leq 3,\; k\in \\{0,1\\}\right\\}.
\]

```python
sage: A = hi.ShiArrangement("A2")
sage: A
Arrangement of 6 hyperplanes of dimension 3 and rank 2
```

The topological zeta function is 

```python
sage: Z = hi.TopologicalZetaFunction(A)
sage: Z
-9/(s + 1) - 3*(s - 2)/((3*s + 2)*(s + 1)) + 3/(s + 1)^2 + 4
```

which is equivalent to
\[
    Z_{\mathcal{S}\mathsf{A}_2}^{\mathrm{top}}(s) = \dfrac{2 - 5s + 2s^2 + 12s^3}{(s+1)^2(3s+2)} .
\]
