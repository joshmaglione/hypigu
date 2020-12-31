# Generating Functions

The main purpose of this package is to explicitly compute the flag Hilbert&ndash;Poincar&#233; series and its specializations like Igusa'a local zeta function. We keep variable letters consistent with Maglione&ndash;Voll; the exceptions is we replace $q^{-s_x}$ with $t_x$, for some label $x$. 

## AnalyticZetaFunction

**Input**:

- a hyperplane arrangement $\mathcal{A}$.

**Output**:

- the analytic zeta function associated to $\mathcal{A}$. 

Given a suitable $\mathfrak{o}$-representation of a $d$-dimensional hyperplane arrangement $\mathcal{A}$, where $\mathfrak{o}$ is a compact discrete valuation ring, the analytic zeta function is defined to be the integral:

\[
    \zeta_{\mathcal{A}(\mathfrak{o})}(\textbf{s}) = \int_{\mathfrak{o}^d} \prod_{x\in\widetilde{\mathcal{L}}(\mathcal{A})} \\| \mathcal{A}_x\\|^{s_x} \\, d\mu. 
\]

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

- a hyperplane arrangement $\mathcal{A}$.

**Output**:

- the atom zeta function associated to $\mathcal{A}$. 

Given a suitable $\mathfrak{o}$-representation of a $d$-dimensional hyperplane arrangement $\mathcal{A}$, where $\mathfrak{o}$ is a compact discrete valuation ring, the atom zeta function is defined to be the integral:

\[
    \zeta_{\mathcal{A}(\mathfrak{o})}^{\mathrm{at}}(\textbf{s}) = \int_{\mathfrak{o}^d} \prod_{L\in\mathcal{A}(\mathfrak{o})} |L|^{s_L} \\, d\mu. 
\]

## CombinatorialSkeleton

**Input**:

- a hyperplane arrangement $\mathcal{A}$.

**Output**:

- the combinatorial skeleton associated to $\mathcal{A}$. 

The combinatorial skeleton of $\mathcal{A}$ is defined to be:

\[
    HP_{\mathcal{A}}^{\mathrm{sk}} (Y, T) 
    = \sum_{F\in\Delta(\widetilde{\mathcal{L}}(\mathcal{A}))} \pi_F(Y) \left(\dfrac{T}{1 - T}\right)^{|F|}.
\]

## FlagHilbertPoincareSeries

**Input**:

- a hyperplane arrangement $\mathcal{A}$.

**Output**:

- the flag Hilbert&ndash;Poincar&#233; series associated to $\mathcal{A}$. 

The flag Hilbert&ndash;Poincar&#233; series of $\mathcal{A}$ is defined to be:

\[
    HP_{\mathcal{A}} (Y, \mathbf{T}) 
    = \sum_{F\in\Delta(\widetilde{\mathcal{L}}(\mathcal{A}))} \pi_F(Y) \prod_{x\in F} \frac{T_x}{1 - T_x}.
\]

## IgusaZetaFunction

**Input**:

- a hyperplane arrangement $\mathcal{A}$ or a polynomial $f$.

**Output**:

- Igusa's local zeta function associated to either $\mathcal{A}$ or $f$. 

If a polynomial, $f$, is given, we require that $f$ be the product of linear factors. Symbolic expressions and strings are fine as well, provided SageMath interprets them as a polynomial.

For a compact discrete valuation ring $\mathfrak{o}$ and a polynomial $f\in \mathfrak{o}[X_1,\dots, X_d]$, Igusa's local zeta function associated with $f$ is 

\[
    Z_f(s) = \int_{\mathfrak{o}^d} |f(\mathbf{X})|^s\, d\mu.
\]

If $Q_\mathcal{A}$ is the defining polynomial of a hyperplane arrangement, then Igusa's local zeta function associated $\mathcal{A}$ is $Z_{Q_\mathcal{A}}(s)$. 


