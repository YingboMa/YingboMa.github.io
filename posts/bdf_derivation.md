+++
title = "A Complete Derivation of Quasi-constant Step Size BDF Method"
hascode = true
date = Date(2021, 4, 23)
rss = "We will do a complete derivation of all the formulas necessary for
implementing a complete quasi-constant step size backward differentiation
formulas method."
+++

@def tags = ["math", "differential equations"]

Backward differentiation formula (BDF) method is the most widely used method
to solve stiff ordinary differential equations (ODEs) and differential-algebraic
equations (DAEs). Famous solvers like `ode15s`, `ode15i`, `LSODA`, `CVODE`,
`IDA`, and `DASSL` are all BDFs. Unfortunately, all BDF derivations online are
not complete in the sense that one still needs to do hand derivation for some
formulas to be able to write a robust and production-ready variable step size
BDF implementation.

Recall the definitions of ODEs $u' = f(u, p, t)$ and DAEs $F(u', u, p, t) = 0$.
The idea behind BDFs is to use a polynomial $p(t)$ to interpolate $u(t)$ using
the points $(t_{n+1}, u_{n+1}), (t_{n}, u_{n}), (t_{n-1}, u_{n-1}), ...$, and
solve for the next step $u_{n+1}$ that makes $p'(t+h) = u'_{n+1}$, where $h$ is
the step size.

\section{Constant Step Size BDF}
Let's first solve the simple problem: deriving BDFs assuming the step size $h$
is constant. An order $s$ BDF interpolates $s+1$ points to form a degree $s$
polynomial $p^{(n+1)}_s(t)$ when computing the approximation $u_{n+1}$. Since we
want to compute the next step, we can center the polynomial at the current time
$t_{n}$ and parametrize $t$ by the step size $h$. Thus, we define $q^{(n+1)}_{s}
(c) = p_{s}(t) = p_{s}(t_{n}+ch)$ where $c\in\R$ is the new independent variable
after the change of variable. Also note that we have

$$t- t_{n-i} = (t_{n} + ch) - (t_{n} - ih) = (c+i)h.$$

We can then explicitly construct the $s$-th order interpolant using the [Newton
polynomial interpolation](http://fourier.eng.hmc.edu/e176/lectures/ch7/node4.html):
\begin{align}
&q^{(n+1)}_{s}(c) = p^{(n+1)}_{s}(t) = p^{(n+1)}_{s}(t_{n}+ch) \\
= \;& u_{n+1} + [u_{n+1},u_{n}](t-t_{n+1}) + ... + [u_{n+1},...,u_{n+1-s}](t-t_{n+1})\cdots (t-t_{n+2-s}) \\
= \;& u_{n+1} + \sum_{j=1}^s [u_{n+1},...,u_{n+1-j}](t-t_{n+1})\cdots(t-t_{n-(j-2)}) \\
= \;& u_{n+1} + \sum_{j=1}^s \frac{(c-1) c \cdots (c+j-2)}{j!} h^j j![u_{n+1},...,u_{n+1-j}] \\
= \;& u_{n+1} + \sum_{j=1}^s \frac{1}{j!} \left( \prod_{i=1}^j c+i-2 \right) \nabla_{h}^{(j)} u_{n+1},
\end{align}
where $[...]$ is the Newton's divided difference, and $\nabla$ is the [backward
differentiation operator](https://en.wikipedia.org/wiki/Finite_difference#Higher-order_differences).
We used the relation $h^j j![u_{n+1},...,u_{n+1-j}] = \nabla_{h}^{(j)} u_{n+1}$
in the last step.

Now, we want to set $p'(t+h) = u'_{n+1}$. Note that we have $q'(1) = h p'(t_{n}
+h) = h u'_{n+1}$. Therefore, we need to solve
$$
\sum_{j=1}^s \delta_{j} \nabla^{(j)} u_{n+1} = h u'_{n+1},
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
\sum_{j=1}^s \frac{1}{j} \nabla^{(j)} u_{n+1} = h u'_{n+1}.
$$

\section{Varying the Step Size by Transplanting to Equidistant Grid}

A simple strategy to vary the step size is to simply evaluate the old polynomial
$p$ from the at the equidistant grid $t_{n} - i\tilde{h}$ for $i=0, 1, 2,..., s$
to obtain an interpolated "constant step size history", when changing to a new
step size $\tilde{h}$. During the step changing process, we do not need to
compute $u_{n+1}$. Hence, we only need the polynomial $p^{(n)}(t_{n}+ch)$. This
strategy is simple because evaluating and re-interpolating a polynomial
interpolant are linear transformations of the history, so we can solve for
matrices that performs such transformations instead of working with polynomials
explicitly.

Note that the only step size dependent term is $\nabla_{h}^{(j)} u_{n}$. Let's
define $r = \tilde{h}/h$,
\begin{align}
D = [\nabla_{h} u_{n}, \nabla_{h}^2 u_{n}, ..., \nabla_{h}^s u_{n}] \in \R^{m \times s},
\quad \text{and} \quad \tilde{D} = [\nabla_{\tilde{h}} u_{n}, \nabla_{\tilde{h}}^2 u_{n}, ..., \nabla_{\tilde{h}}^s u_{n}]  \in \R^{m \times s}.
\end{align}
Constructing the polynomial $p^{(n)}_{s}(t_{n}+ch)$ is simple as we only need to
subtract all the indeices in $p^{(n+1)}_{s}(t_{n}+ch)$ by 1, which results in
\begin{align}
p^{(n)}_{s}(t_{n}+ch) = \;& q^{(n)}_{s}(c) = u_{n} + \sum_{j=1}^s [u_{n},...,u_{n-j}](t-t_{n})\cdots(t-t_{n-(j-1)}) \\
= \;& u_{n} + \sum_{j=1}^s \frac{1}{j!} \left( \prod_{i=1}^j c+i-1 \right) \nabla_{h}^{(j)} u_{n}\\
= \;& u_{n} + \sum_{j=1}^s \frac{1}{j!} \left( \prod_{i=1}^j c+i-1 \right) D_{j}\\
= \;& u_{n} + \sum_{j=1}^s \frac{1}{j!} \left( \prod_{i=0}^{j-1} c+i \right) D_{j}.
\end{align}
The interpolating polynomial with the new step size $\tilde{h} = rh$ is then
\begin{align}
\tilde{p}^{(n)}_{s}(t_{n}+\tilde{h}) = \tilde{q}^{(n)}_{s}(c) = u_{n} + \sum_{j=1}^s \frac{1}{j!} \left( \prod_{i=0}^{j-1} c+i \right) \tilde{D}_{j}.
\end{align}
With the interpolation condition, we force $p^{(n)}_{s}(t_{n}+c\tilde{h}) =
p^{(n)}_{s}(t_{n}+crh) = \tilde{p}^{(n)}_{s}(t_{n}+c\tilde{h})$ for $c = 0, -1,
..., -s$. Note that when $c=0$, the interpolation condition simplifies to
$\tilde{p}^{(n)}_{s}(t_{n}) = u_{n} = p^{(n)}_{s}(t_{n})$ which always holds.
Therefore, we only need to solve
\begin{align}
\sum_{j=1}^s D_{i,j} \frac{1}{j!} \left( \prod_{i=0}^{j-1} cr+i \right) =
\sum_{j=1}^s \tilde{D}_{i,j} \frac{1}{j!} \left( \prod_{i=0}^{j-1} c+i \right)
\end{align}
for $\tilde{D}$ with $c = -1, -2, .., -s$ and $i = 1, ..., m$. When express the
above system of equations in terms of matrix, there is
\begin{align}
D R = \tilde{D} U,
\end{align}
where
\begin{align}
R_{jk} = \frac{1}{j!} \left( \prod_{i=0}^{j-1} i-k r \right), \quad\text{and}
\quad U_{jk} = \frac{1}{j!} \left( \prod_{i=0}^{j-1} i-k \right).
\end{align}
Note that $k = -c$ as $k=1,2,...,s$. We end this section by remarking that
$R^2=I$, so
\begin{align}
\tilde{D} = D (R U).
\end{align}
