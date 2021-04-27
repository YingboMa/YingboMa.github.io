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
polynomial $p_s(t)$. Since we want to compute the next step, we can center the
polynomial at the current time $t_{n}$ and parametrize $t$ by the step size $h$.
Thus, we define $q_{s}(c) = p_{s}(t) = p_{s}(t_{n}+ch)$ where $c\in\R$ is the
new independent variable after the change of variable. Also note that we have

$$t- t_{n-i} = (t_{n} + ch) - (t_{n} - ih) = (c+i)h.$$

We can then explicitly construct the $s$-th order interpolant using the [Newton
polynomial interpolation](http://fourier.eng.hmc.edu/e176/lectures/ch7/node4.html):
\begin{align}
&q_{s}(c) = p_{s}(t) = p_{s}(t_{n}+ch) \\
= \;& [u_{n+1}] + [u_{n+1},u_{n}](t-t_{n+1}) + ... + [u_{n+1},...,u_{n+1-s}](t-t_{n+1})\cdots (t-t_{n+2-s}) \\
= \;& [u_{n+1}] + \sum_{j=1}^s [u_{n+1},...,u_{n+1-j}](t-t_{n+1})\cdots(t-t_{n-(j-2)}) \\
= \;& [u_{n+1}] + \sum_{j=1}^s \frac{(c-1) c \cdots (c+j-2)}{j!} h^j j![u_{n+1},...,u_{n+1-j}] \\
= \;& [u_{n+1}] + \sum_{j=1}^s \frac{1}{j!} \left( \prod_{i=1}^j c+i-2 \right) \nabla_{h}^{(j)} u_{n+1},
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
$p$ at the equidistant grid $t_{n} - i\tilde{h}$ for $i=1, 2,..., s-2$ to
obtain an interpolated "constant step size history", when changing to a new step
size $\tilde{h}$. This strategy is simple because evaluating and
re-interpolating a polynomial interpolant are linear transformations of the
history, so we can solve for matrices that performs such transformations instead
of working with polynomials explicitly.

Recall that $q_{s}(c) = p_{s}(t) = p_{s}(t_{n}+ch)$, we will work on the
transformed $c$-coordinate for simplicity. Let's denote $\tilde{q}_{s}$ as the
re-interpolated polynomial with the new step size $\tilde{h}$. Note that we have
\begin{align}
q_{s}(c) = [u_{n+1}] + \sum_{j=1}^s \frac{1}{j!} \left( \prod_{i=1}^j c+i-2 \right) \nabla_{h}^{(j)} u_{n+1}, \\
\tilde{q}_{s}(c) = [u_{n+1}] + \sum_{j=1}^s \frac{1}{j!} \left( \prod_{i=1}^j c+i-2 \right) \nabla_{\tilde{h}}^{(j)} u_{n+1}.
\end{align}
The only step size dependent term is $\nabla_{h}^{(j)} u_{n+1}$. Let's define
\begin{align}
D = [\nabla_{h} u_{n+1}, \nabla_{h}^2 u_{n+1}, ...,  , \nabla_{h}^s u_{n+1}]
\quad \text{and} \quad \tilde{D} = [\nabla_{\tilde{h}} u_{n+1}, \nabla_{\tilde{h}}^2 u_{n+1}, ...,  , \nabla_{\tilde{h}}^s u_{n+1}].
\end{align}
With the interpolation condition, we force $q_{s}(c) = \tilde{q}_{s}(c)$ for $c
= -1, 0, ..., s-2$:
\begin{align}
\sum_{j=1}^s \frac{1}{j!} \left( \prod_{i=1}^j c+i-2 \right) \nabla_{h}^{(j)} u_{n+1} =
\sum_{j=1}^s \frac{1}{j!} \left( \prod_{i=1}^j c+i-2 \right) \nabla_{\tilde{h}}^{(j)} u_{n+1}.
\end{align}
