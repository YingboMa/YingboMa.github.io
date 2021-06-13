+++
title = "A Complete Derivation of Quasi-constant Step Size BDF Method"
hascode = true
date = Date(2021, 4, 23)
rss = "We will do a complete derivation of all the formulas necessary for
implementing a complete quasi-constant step size backward differentiation
formulas method."
+++

@def tags = ["math", "differential equations"]

The backward differentiation formula (BDF) method is the most widely used method
to solve stiff ordinary differential equations (ODEs) and differential-algebraic
equations (DAEs). Famous stiff solvers like `ode15s`, `ode15i`, `LSODE`, `LSODA`,
`CVODE`, `IDA`, and `DASSL` all use the BDF method. Unfortunately, all BDF
derivations online are not complete, in the sense that one still needs to do
hand derivation for some formulas to fully understand the algorithm and be able
to write a robust and production-ready variable step size BDF implementation.
This blog post aims to bridge this gap by providing a complete derivation of all
the formulas needed for a production-ready quasi-constant step size BDF method.

Recall the definitions of ODEs $u' = f(u, p, t)$ and DAEs $F(u', u, p, t) = 0$.
The idea behind BDFs is to use a polynomial $p(t)$ to interpolate $u(t)$ using
the points $(t_{n+1}, u_{n+1}), (t_{n}, u_{n}), (t_{n-1}, u_{n-1}), ...$, and
solve for the next step $u_{n+1}$ that makes $p'(t+h) = u'_{n+1}$, where $h$ is
the step size.

\section{Constant step size BDF}
Let's first solve the simple problem: deriving BDFs assuming the step size $h$
is constant. An order $s$ BDF interpolates $s+1$ points to form a degree $s$
polynomial $p_{s,n+1}(t)$ when computing the approximation $u_{n+1}$. Since we
want to compute the next step, we can center the polynomial at the current time
$t_{n}$ and parametrize $t$ by the step size $h$. Thus, we define $q_{s,n+1}
(c) = p_{s,n+1}(t) = p_{s,n+1}(t_{n}+ch)$ where $c\in\R$ is the new independent
variable after the change of variable. Also note that we have

$$t- t_{n-i} = (t_{n} + ch) - (t_{n} - ih) = (c+i)h.$$

We can then explicitly construct the $s$-th order interpolant using the [Newton
polynomial interpolation](http://fourier.eng.hmc.edu/e176/lectures/ch7/node4.html):
\begin{align}
&q_{s,n+1}(c) = p_{s,n+1}(t) = p_{s,n+1}(t_{n}+ch) \\
= \;& u_{n+1} + [u_{n+1},u_{n}](t-t_{n+1}) + ... + [u_{n+1},...,u_{n+1-s}](t-t_{n+1})\cdots (t-t_{n+2-s}) \\
= \;& u_{n+1} + \sum_{j=1}^s [u_{n+1},...,u_{n+1-j}](t-t_{n+1})\cdots(t-t_{n-(j-2)}) \\
= \;& u_{n+1} + \sum_{j=1}^s \frac{(c-1) c \cdots (c+j-2)}{j!} h^j j![u_{n+1},...,u_{n+1-j}] \\
= \;& u_{n+1} + \sum_{j=1}^s \frac{1}{j!} \left( \prod_{i=1}^j c+i-2 \right) \nabla_{h}^{j} u_{n+1},
\end{align}
where $[...]$ denotes Newton's divided difference, and $\nabla$ is the [backward
differentiation operator](https://en.wikipedia.org/wiki/Finite_difference#Higher-order_differences).
We used the relation $h^j j![u_{n+1},...,u_{n+1-j}] = \nabla_{h}^{j} u_{n+1}$
in the last step.

Now, we want to set $p'(t+h) = u'_{n+1}$. Note that we have $q'(1) = h p'(t_{n}
+h) = h u'_{n+1}$. Therefore, we need to solve
$$
\sum_{j=1}^s \delta_{j} \nabla^{j} u_{n+1} = h u'_{n+1},
$$
where
\begin{align}
\delta_{j} &= \frac{1}{j!} \dd{}{c} \prod_{i=1}^j c-i-2 \Big|_{c=1} \\
&= \frac{1}{j!} \sum_{i=1}^j \prod_{j=1 \text{ and } j\ne i}^j c + i - 2 \Big|_{c=1} \\
&= \frac{1}{j!} \sum_{i=1}^j \prod_{j=1 \text{ and } j\ne i}^j j-1 \quad \text{note that the product vanishes when $i\ne 1$}\\
&= \frac{1}{j!} \prod_{j=2}^j j-1 = \frac{1}{j!} (j-1)! = \frac{1}{j}.
\end{align}

Together, the constant step size BDF is
$$
\sum_{j=1}^s \frac{1}{j} \nabla^{j} u_{n+1} = h u'_{n+1}.
$$

\section{Varying the step size by transplanting to equidistant grid}

A simple strategy to vary the step size is to simply evaluate the old polynomial
$p$ from the at the equidistant grid $t_{n} - i\tilde{h}$ for $i=0, 1, 2,..., s$
to obtain an interpolated "constant step size history", when changing to a new
step size $\tilde{h}$. During the step changing process, we do not need to
compute $u_{n+1}$. Hence, we only need the polynomial $p_{s,n}(t_{n}+ch)$. This
strategy is simple because evaluating and re-interpolating a polynomial
interpolant are linear transformations of the history, so we can solve for
matrices that performs such transformations instead of working with polynomials
explicitly.

Note that the only step size dependent term is $\nabla_{h}^{j} u_{n}$. Let's
define $r = \tilde{h}/h$,
\begin{align}
D = [\nabla_{h} u_{n}, \nabla_{h}^2 u_{n}, ..., \nabla_{h}^s u_{n}] \in \R^{m \times s},
\quad \text{and} \quad \tilde{D} = [\nabla_{\tilde{h}} u_{n}, \nabla_{\tilde{h}}^2 u_{n}, ..., \nabla_{\tilde{h}}^s u_{n}]  \in \R^{m \times s}.
\end{align}
Constructing the polynomial $p_{s,n}(t_{n}+ch)$ is simple as we only need to
subtract the indices in $p_{s,n+1}(t_{n}+ch)$ by 1 appropriately:
\begin{align}
p_{s,n}(t_{n}+ch) = q_{s,n}(c) = \;& u_{n} + \sum_{j=1}^s [u_{n},...,u_{n-j}](t-t_{n})\cdots(t-t_{n-(j-1)}) \\
= \;& u_{n} + \sum_{j=1}^s \frac{1}{j!} \left( \prod_{i=1}^j c+i-1 \right) \nabla_{h}^{j} u_{n}\\
= \;& u_{n} + \sum_{j=1}^s \frac{1}{j!} \left( \prod_{i=1}^j c+i-1 \right) D_{j}\\
= \;& u_{n} + \sum_{j=1}^s \frac{1}{j!} \left( \prod_{i=0}^{j-1} c+i \right) D_{j}.
\end{align}
The interpolating polynomial with the new step size $\tilde{h} = rh$ is then
\begin{align}
\tilde{p}_{s,n}(t_{n}+\tilde{h}) = \tilde{q}_{s,n}(c) = u_{n} + \sum_{j=1}^s \frac{1}{j!} \left( \prod_{i=0}^{j-1} c+i \right) \tilde{D}_{j}.
\end{align}
With the interpolation condition, we force $p_{s,n}(t_{n}+c\tilde{h}) =
p_{s,n}(t_{n}+crh) = \tilde{p}_{s,n}(t_{n}+c\tilde{h})$ for $c = 0, -1,
..., -s$. Note that when $c=0$, the interpolation condition simplifies to
$\tilde{p}_{s,n}(t_{n}) = u_{n} = p_{s,n}(t_{n})$ which always holds.
Therefore, we only need to solve
\begin{align}
\sum_{j=1}^s D_{i,j} \frac{1}{j!} \left( \prod_{i=0}^{j-1} cr+i \right) =
\sum_{j=1}^s \tilde{D}_{i,j} \frac{1}{j!} \left( \prod_{i=0}^{j-1} c+i \right)
\end{align}
for $\tilde{D}$ with $c = -1, -2, .., -s$ and $i = 1, ..., m$. The system of
equations can be expressed in terms of matrices:
\begin{align}
D R = \tilde{D} U,
\end{align}
where
\begin{align}
R_{jk} = \frac{1}{j!} \left( \prod_{i=0}^{j-1} i-k r \right), \quad\text{and}
\quad U_{jk} = \frac{1}{j!} \left( \prod_{i=0}^{j-1} i-k \right).
\end{align}
Note that $c = -k$ as $k=1,2,...,s$. We end this section by remarking that
$U^2=I$.
@@colbox-blue
```julia
julia> using LinearAlgebra

julia> U(s) = [prod(i->i-k, 0:j-1)//factorial(j) for j in 1:s, k in 1:s]
U (generic function with 1 method)

julia> U(5)^2 == I
true
```
@@
Therefore, the divided differences with $\tilde{h}$ is simply
\begin{align}
\tilde{D} = D (R U).
\end{align}

\section{Updating backward differences}
A natural choice for the predictor is then $u^{(0)}_{n+1} = p_{s,n}(t_{n}+h) = u_{n} + \sum_{j=1}
^s D_{j}$. From the definition of the backward difference, we have
\begin{align}
\nabla_{h}^{s+1} u_{n+1} &= \nabla_{h}^{s} u_{n+1} - \nabla_{h}^{s} u_{n}\\
&= \nabla_{h}^{s-1} u_{n+1} - \nabla_{h}^{s-1} u_{n} - \nabla_{h}^{s} u_{n} \\
&= \nabla_{h}^{s-2} u_{n+1} - \nabla_{h}^{s-2} u_{n} - \nabla_{h}^{s-1} u_{n} - \nabla_{h}^{s} u_{n} \\
&= u_{n+1} - u_{n} - ... - \nabla_{h}^{s-2} u_{n} - \nabla_{h}^{s-1} u_{n} - \nabla_{h}^{s} u_{n} \\
&= u_{n+1} - u^{(0)}_{n+1}.
\end{align}
Therefore, the $s+1$-th order backward difference of the next step is
simply the difference between the corrector $u_{n+1}$ and the predictor $u^{(0)}
_{n+1}$. Also, note that from
\begin{align}
\nabla_{h}^{j} u_{n+1} = \nabla_{h}^{j+1} u_{n+1} + \nabla_{h}^{j} u_{n},
\end{align}
where $j\in\mathbb{N}$, we can compute all the lower order backward differences
of the next step, and
\begin{align}
\nabla_{h}^{s+2} u_{n+1} = \nabla_{h}^{s+1} u_{n+1} - \nabla_{h}^{s+1} u_{n}.
\end{align}
We can use the above updating relations to simplify the backward differences:
\begin{align}
\nabla^s u_{n+1} &= \nabla^{s+1} u_{n+1} + \nabla^{s} u_{n} = (u_{n+1}-u_{n+1}^{(0)}) + \nabla^{s} u_{n} \\
\nabla^{s-1} u_{n+1} &= \nabla^{s} u_{n+1} + \nabla^{s-1} u_{n} = (u_{n+1}-u_{n+1}^{(0)}) + \nabla^{s} u_{n}  + \nabla^{s-1} u_{n} \\
&\vdots \\
\nabla^{j} u_{n+1} &= \nabla^{j+1} u_{n+1} + \nabla^{j-1} u_{n} = (u_{n+1}-u_{n+1}^{(0)}) + \sum_{k=j}^{s}\nabla^{k} u_{n}.
\end{align}
Therefore, BDF becomes
\begin{align}
\sum_{j=1}^s \frac{1}{j} \nabla^{j} u_{n+1} &= \sum_{j=1}^s \left[\frac{1}{j}\left(
(u_{n+1}-u_{n+1}^{(0)}) + \sum_{k=j}^{s}\nabla^{k} u_{n}\right)\right] \\
&=\gamma_{s} (u_{n+1}-u_{n+1}^{(0)}) + \sum_{j=1}^s \frac{1}{j}\sum_{k=j}^{s} \nabla^{k} u_{n},
\end{align}
where $\gamma_{j} = \sum_{k=1}^j \frac{1}{j}$. This naïve approach contains a
lot of redundant computation, and the computational complexity is $O(s^2 m)$ for
$u_{n}\in \R^m$. We can do a lot better if reorder the summation:
\begin{align}
&\sum_{j=1}^s \frac{1}{j}\sum_{k=j}^{s} \nabla^{k} u_{n}
= \sum_{j=1}^s\sum_{k=j}^{s} \frac{1}{j} \nabla^{k} u_{n}
= \sum_{j=1}^s\sum_{k=1 \;\land\; k\ge j}^{s} \frac{1}{j} \nabla^{k} u_{n}\\
=&\sum_{k=1}^s\sum_{j=1 \;\land\; j\le k}^{s} \frac{1}{j} \nabla^{k} u_{n}
\quad\text{note the interchange of summation}\\
=&\sum_{k=1}^s\sum_{j=1}^{k} \frac{1}{j} \nabla^{k} u_{n}
= \sum_{k=1}^s\left(\sum_{j=1}^{k} \frac{1}{j}\right) \nabla^{k} u_{n}
= \sum_{k=1}^s \gamma_{k} \nabla^{k} u_{n}.
\end{align}
This expression can be evaluated in $O(sm)$ time.

Putting everything together, the simplified BDF is
\begin{align}
\gamma_{s} (u_{n+1}-u_{n+1}^{(0)}) +
\sum_{k=1}^s \gamma_{k} \nabla^{k} u_{n} = h u'_{n+1},
\end{align}
and we can use a nonlinear equation solver to solve for $u_{n+1}$. We will omit
the details of writing a stable and efficient nonlinear solver here as they are
independent from the BDF method itself.

\section{Local truncation error}
By the standard result on polynomial interpolation, we know that
\begin{align}
u(t) - p_{s,n+1}(t) = \frac{u^{(s+1)}(\xi_{t})}{(s+1)!} w(t),
\end{align}
where $w(t) = \prod_{j=-1}^{s-1}(t - t_{n-j})$. We assume that all the history
is completely accurate and compute the error of $p'_{s,n+1}(t_{n+1})$:
\begin{align}
u'(t_{n+1}) = p'_{s,n+1}(t_{n+1}) + \frac{u^{(s+1)}(\xi_{t})}{(s+1)!} w'(t_{n+1}).
\end{align}
Note that we don't have the term with $w(t_{n+1})$ because it goes to $0$ by
the construction of $w$. Expanding $w'(t_{n+1})$:
\begin{align}
w'(t_{n+1}) &= \sum_{j=-1}^{s-1} \left[(t_{n+1}-t_{n-j})' \prod_{k=-1\;\land\; k\ne j}^{s-1} (t_{n+1}-t_{n-k}) \right] \\
&= \sum_{j=-1}^{s-1} \prod_{k=-1\;\land\; k\ne j}^{s-1} (t_{n+1}-t_{n-k}) \quad
\text{the product only doesn't vanish when $j=-1$} \\
&= \prod_{k=0}^{s-1} (t_{n+1}-t_{n-k}).
\end{align}
Together, BDF with the leading error term is then
\begin{align}
hu'(t_{n+1}) &= hp'_{s,n+1}(t_{n+1}) + h\frac{u^{(s+1)}(\xi_{t})}{(s+1)!} \prod_{k=0}^{s-1} (t_{n+1}-t_{n-k}) \\
&= q'_{s,n+1}(t_{n+1}) + h\frac{u^{(s+1)}(\xi_{t})}{(s+1)!} (h\cdot 2h
\cdot 3h ... sh) \\
&= q'_{s,n+1}(t_{n+1}) + h\frac{u^{(s+1)}(\xi_{t})}{(s+1)!} s!h^s \\
&= q'_{s,n+1}(t_{n+1}) + \frac{u^{(s+1)}(\xi_{t})}{s+1} h^{s+1} \\
&= \sum_{j=1}^s \frac{1}{j} \nabla^{j} u_{n+1} + \frac{1}{s+1}\nabla_{h}^{s+1} u_{n+1} + O(h^{s+2}).
\end{align}
Hence, the error estimate for the $s$-th order BDF is
\begin{align}
\frac{1}{s+1}\nabla_{h}^{s+1} u_{n+1},
\end{align}
and this quantity can easily be computed from the updated backward differences.
With an appropriate controller for the step size and order, a quasi-constant
step size BDF method can be implemented in its totality.
