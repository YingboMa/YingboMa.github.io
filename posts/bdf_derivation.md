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
the historical points $(t_{n}, u_{n}), (t_{n-1}, u_{n-1}), ...$, and solve for
next step $u_{n+1}$ that makes $p'(t+h) = u'_{n+1}$ where $h$ is the step size.

\section{Constant Step Size BDF}
Let's first solve the simple problem: deriving BDFs assuming the step size $h$
is constant. An order $s$ BDF interpolates $s+1$ points to form a degree $s$
polynomial $p_s(t)$. Since we want to compute the next step, we can center the
polynomial at the current time $t_{n}$ and parametrize $t$ by the step size $h$.
Thus, we define $q_{s}(c) = p_{s}(t) = p_{s}(t_{n}+ch)$ where $c\in\R$ is the
new independent variable after the change of variable. Also note that we have

$$t- t_{n-i} = (t_{n} + ch) - (t_{n} - ih) = (c+i)h.$$

We can then explicitly construct the interpolant using the Newton polynomial
interpolation:
\begin{align}
&q_{s}(c) = p_{s}(t) = p_{s}(t_{n}+ch) \\
= \;& [u_{n+1}] + [u_{n+1},u_{n}](t-t_{n+1}) + ... + [u_{n+1},...,u_{n+1-s}](t-t_{n+1})\cdots (t-t_{n+2-s}) \\
= \;& [u_{n+1}] + \sum_{j=1}^s [u_{n+1},...,u_{n+1-j}](t-t_{n+1})\cdots(t-t_{n-(j-2)}) \\
= \;& [u_{n+1}] + \sum_{j=1}^s \frac{(c-1) c \cdots (c+j-2)}{j!} h^j j![u_{n+1},...,u_{n+1-j}] \\
= \;& [u_{n+1}] + \sum_{j=1}^s \frac{1}{j!} \left( \prod_{i=1}^j c+i-2 \right) \Delta^{(j)} u_{n+1},
\end{align}
where $[...]$ is the Newton's divided difference, and $\Delta$ is the backward
differentiation operator. We used the relation $h^j j![u_{n+1},...,u_{n+1-j}] =
\Delta^{(j)} u_{n+1}$ in the last step.

Now, we want to set $p'(t+h) = u'_{n+1}$. Note that we have $q'(1) = h p'(t_{n}
+h) = h u'_{n+1}$. Therefore, we need to solve
$$
\sum_{j=1}^s \delta_{j} \Delta^{(j)} u_{n+1} = h u'_{n+1},
$$
where
\begin{align}
\delta_{j} &= \frac{1}{j!} \dd{}{c} \prod_{i=1}^j c-i-2 \Big|_{c=1} \\
&= \frac{1}{j!} \sum_{i=1}^j \prod_{j=1 \text{ and } j\ne i}^j c + i - 2 \Big|_{c=1} \\
&= \frac{1}{j!} \sum_{i=1}^j \prod_{j=1 \text{ and } j\ne i}^j j-1 \quad \text{note that the product vanishes when $i\ne 1$}\\
&= \frac{1}{j!} \prod_{j=2}^j j-1 = \frac{1}{j!} (j-1)! = \frac{1}{j}.
\end{align}

Hence, the constant step size BDF is
$$
\sum_{j=1}^s \frac{1}{j} \Delta^{(j)} u_{n+1} = h u'_{n+1}.
$$
