+++
title = "Taylor series and trees"
hascode = true
date = Date(2025, 9, 19)
rss = "Taylor series and trees"
+++

@def tags = ["math"]

\section{Runge–Kutta, Taylor truncation, and “tree” bookkeeping in elementary index notation}

We build everything from two ordinary calculus rules and nothing else:
- Chain rule: for smooth $g: \mathbb{R}^d\to\mathbb{R}$, $\tfrac{\mathrm{d}}{\mathrm{d}t} g(u(t)) = \sum_{b=1}^d \tfrac{\partial g(u(t))}{\partial u_b}\, u_b'(t)$.
- Product rule: $\tfrac{\mathrm{d}}{\mathrm{d}t} \prod_{r=1}^m G_r(t) = \sum_{s=1}^m \Bigl( \prod_{r\ne s} G_r(t) \Bigr) G_s'(t)$.

We apply these two rules repeatedly to the autonomous ODE $u_i'(t) = f_i(u(t))$ to obtain a symbolic Taylor expansion and an indexing set that later becomes “trees.” No advanced machinery, and no presupposed form.

\section{1) Elementary construction: Taylor of $u$ and $f$ $\Rightarrow$ Taylor of $f(u(t_0+h))$ (order $k$)}

We fix an order $k\ge 1$ and work only with Taylor polynomials and coefficient extraction, all modulo $h^{k+1}$. Throughout, $u_0 = u(t_0)$, $U(h) = (U_1(h),\dots,U_d(h))^\top$, and $U^{[q]}\in\mathbb{R}^d$ denotes the vector of $q$-th Taylor coefficients at $t_0$.

\section{1.1 Build the Taylor polynomial of $u(t_0+h)$}

Write the degree-$k$ Taylor polynomial of $u$ at $t_0$ with explicit truncation:
$$
  u_i(t_0+h) \equiv u_{0,i} + \sum_{q=1}^{k} U_i^{[q]} h^q \pmod{h^{k+1}},\qquad U_i^{[q]} = \frac{1}{q!}\, u_i^{(q)}(t_0).
$$
Collecting components,
$$
  U(h) \;=\; u_0 + \sum_{q=1}^{k} U^{[q]} h^q \in \mathbb{R}^d[h]/(h^{k+1}).
$$
Notation: $U_a(h)$ is the $a$-th component of $U(h)$; $u_{0,a}$ is the $a$-th component of $u_0$.

\section{1.2 Build the Taylor polynomial of $f$ at $u_0$}

For $f: \mathbb{R}^d\to\mathbb{R}^d$ with $C^{k}$ smoothness, the degree-$k$ Taylor polynomial around $u_0$ is the finite polynomial (componentwise)
\begin{equation}\label{eq:taylor-f}
  f_i(u_0 + y) \equiv \sum_{m=0}^{k} \frac{1}{m!}
  \sum_{a_1=1}^d\!\cdots\!\sum_{a_m=1}^d
  \frac{\partial^{m} f_i(u_0)}{\partial u_{a_1}\cdots\partial u_{a_m}}\, y_{a_1} \cdots y_{a_m} \pmod{\lVert y\rVert^{k+1}}.
\end{equation}
Construction from the scalar Taylor expansion: fix the increment $y\in\mathbb{R}^d$ and define the scalar function $\phi_i(t):=f_i(u_0 + t\,y)$. Apply the 1D Taylor expansion of $\phi_i$ at $t=0$ up to order $k$:
\begin{equation}\label{eq:phi-taylor}
  \phi_i(t) \equiv \sum_{m=0}^{k} \frac{\phi_i^{(m)}(0)}{m!}\, t^{m} \pmod{t^{k+1}}.
\end{equation}
Compute the first derivatives at $t=0$ by direct differentiation of $\phi_i(t)$:
$$
  \phi_i'(0) = \sum_{a} \frac{\partial f_i(u_0)}{\partial u_a}\, y_a,\qquad
  \phi_i''(0) = \sum_{a,b} \frac{\partial^2 f_i(u_0)}{\partial u_a\partial u_b}\, y_a y_b,\qquad
  \phi_i^{(3)}(0) = \sum_{a,b,c} \frac{\partial^3 f_i(u_0)}{\partial u_a\partial u_b\partial u_c}\, y_a y_b y_c.
$$
Continuing in this way, each differentiation introduces one additional $y$ and one additional partial derivative index, yielding for order $m$ the term
$$
  \phi_i^{(m)}(0) = \sum_{a_1,\dots,a_m} \frac{\partial^{m} f_i(u_0)}{\partial u_{a_1}\cdots\partial u_{a_m}}\, y_{a_1}\cdots y_{a_m}.
$$
Insert these into the scalar Taylor polynomial for $\phi_i(t)$ and then evaluate at $t=1$ to obtain exactly the multivariable polynomial above for $f_i(u_0+y)$, now constructed directly from the scalar Taylor expansion along the line $t\mapsto u_0+t y$.

Elementary differentials (we will reference these throughout):
\begin{equation}\label{eq:elem-bullet}
  F_i(\bullet)(u) = f_i(u),
\end{equation}
\begin{equation}\label{eq:elem-branch}
  F_i\bigl([\tau_1,\dots,\tau_m]\bigr)(u) = \sum_{a_1,\dots,a_m}
   \biggl( \frac{\partial^m f_i(u)}{\partial u_{a_1}\cdots\partial u_{a_m}} \biggr)
   \prod_{r=1}^{m} F_{a_r}(\tau_r)(u).
\end{equation}

\section{1.3 Compose and read off coefficients}

Substitute $y = U(h) - u_0$ into the polynomial for $f$ to obtain the degree-$k$ Taylor polynomial of the composition:
\begin{equation}\label{eq:compose}
  f_i\bigl( u(t_0+h) \bigr)
  \equiv \sum_{m=0}^{k} \frac{1}{m!}
  \sum_{a_1=1}^d\!\cdots\!\sum_{a_m=1}^d
  \frac{\partial^{m} f_i(u_0)}{\partial u_{a_1}\cdots\partial u_{a_m}}
  \Bigl( U_{a_1}(h) - u_{0,a_1} \Bigr) \cdots \Bigl( U_{a_m}(h) - u_{0,a_m} \Bigr) \pmod{h^{k+1}}.
\end{equation}
Here the increment is defined once and for all by
\begin{equation}\label{eq:increment}
  y(h) := u(t_0+h) - u_0 \equiv U(h) - u_0 \pmod{h^{k+1}},\qquad y_a(h) = U_a(h) - u_{0,a} = \sum_{q=1}^{k} U_a^{[q]} h^q.
\end{equation}
Thus every occurrence of $U_{a_r}(h) - u_{0,a_r}$ is simply $y_{a_r}(h)$ written componentwise.

Step-by-step coefficient extraction from the product of sums (no skips):
1) Write the $m$-fold product explicitly and distribute one summand from each factor:
\begin{equation}\label{eq:product-expand}
  \prod_{r=1}^{m} \Bigl( \sum_{q'=1}^{k} U_{a_r}^{[q']} h^{q'} \Bigr)
  \,=\, \sum_{q_1=1}^{k}\cdots\sum_{q_m=1}^{k} \Bigl( \prod_{r=1}^{m} U_{a_r}^{[q_r]} \Bigr) h^{\, q_1+\cdots+q_m}.
\end{equation}
Why this identity holds (distributivity, step by step):
- Case $m=2$:
$$
  \Bigl( \sum_{q_1=1}^{k} U_{a_1}^{[q_1]} h^{q_1} \Bigr) \Bigl( \sum_{q_2=1}^{k} U_{a_2}^{[q_2]} h^{q_2} \Bigr)
  \,=\, \sum_{q_1=1}^{k} \sum_{q_2=1}^{k} U_{a_1}^{[q_1]} U_{a_2}^{[q_2]}\, h^{\, q_1+q_2}.
$$
  This is just the usual distributive law: every pair $(q_1,q_2)$ picks one summand from each factor.
- Induction step: suppose for $m-1$ factors we have
$$
  \prod_{r=1}^{m-1} \Bigl( \sum_{q'=1}^{k} U_{a_r}^{[q']} h^{q'} \Bigr)
  \,=\, \sum_{q_1=1}^{k}\!\cdots\!\sum_{q_{m-1}=1}^{k} \Bigl( \prod_{r=1}^{m-1} U_{a_r}^{[q_r]} \Bigr) h^{\, q_1+\cdots+q_{m-1}}.
$$
  Multiply both sides by the $m$-th sum and distribute once more:
$$
  \Bigl( \sum_{q_m=1}^{k} U_{a_m}^{[q_m]} h^{q_m} \Bigr) \sum_{q_1=1}^{k}\!\cdots\!\sum_{q_{m-1}=1}^{k} \Bigl( \prod_{r=1}^{m-1} U_{a_r}^{[q_r]} \Bigr) h^{\, q_1+\cdots+q_{m-1}}
  \,=\, \sum_{q_1=1}^{k}\!\cdots\!\sum_{q_m=1}^{k} \Bigl( \prod_{r=1}^{m} U_{a_r}^{[q_r]} \Bigr) h^{\, q_1+\cdots+q_m}.
$$
  This is exactly the stated formula for $m$ factors.

2) Apply the coefficient extractor to pick out the $h^q$-terms:
\begin{equation}\label{eq:coeff-pick}
  \bigl[ h^q \bigr] \; \prod_{r=1}^{m} \sum_{q'=1}^{k} U_{a_r}^{[q']} h^{q'}
  \,=\, \sum_{\substack{q_1,\dots,q_m\in\{1,\dots,k\}\\ q_1+\cdots+q_m = q}} \prod_{r=1}^{m} U_{a_r}^{[q_r]}.
\end{equation}
3) Since $U_{a_r}^{[q']}$ is zero outside $q'\in\{1,\dots,k\}$ by definition, we can state the constraints simply as $q_r\ge 1$ and $q_1+\cdots+q_m=q$:
\begin{equation}\label{eq:coeff-cauchy}
  \boxed{\;\bigl[ h^q \bigr] \; \prod_{r=1}^{m} \sum_{q'=1}^{k} U_{a_r}^{[q']} h^{q'}
  \,=\, \sum_{\substack{q_1+\cdots+q_m = q \\
                        q_r \ge 1\ \text{for } r=1,\dots,m}}
       \prod_{r=1}^{m} U_{a_r}^{[q_r]}\;}.
\end{equation}
Therefore, collecting coefficients defines (for $q=0,\dots,k$)
\begin{equation}\label{eq:F-series}
  f_i\bigl( u(t_0+h) \bigr) \equiv \sum_{q=0}^{k} F_i^{[q]} h^q,\qquad
  F_i^{[0]} = f_i(u_0).
\end{equation}
\begin{equation}\label{eq:Fq-U}
  F_i^{[q]} = \sum_{m=1}^{q} \frac{1}{m!}
    \sum_{a_1=1}^d\!\cdots\!\sum_{a_m=1}^d
    \frac{\partial^{m} f_i(u_0)}{\partial u_{a_1}\cdots\partial u_{a_m}}
    \sum_{\substack{q_1+\cdots+q_m = q \\
                    q_r \ge 1\ \text{for } r=1,\dots,m}}
    \prod_{r=1}^{m} U_{a_r}^{[q_r]}\,.
\end{equation}
Matching coefficients in the identity $u'(t_0+h) \equiv \sum_{q\ge 0} (q+1)\, U_i^{[q+1]} h^{q} = f_i(u(t_0+h)) \equiv \sum_{q\ge 0} F_i^{[q]} h^{q}$ yields, for all $q\ge 0$,
\begin{equation}\label{eq:U-F}
  F_i^{[q]} \;=\; (q+1)\, U_i^{[q+1]},\qquad\text{equivalently}\qquad
  U_i^{[q]} \;=\; \frac{1}{q}\, F_i^{[q-1]}\ \ (q\ge 1).
\end{equation}

Combining this with the Taylor polynomial of $u$ from §1.1 gives the labeled $u$-series we will reference later:
\begin{equation}\label{eq:u-F-sum-reprise}
  u_i(t_0+h) \equiv u_{0,i} + \sum_{q=1}^{k} \frac{h^q}{q}\, F_i^{[q-1]} \pmod{h^{k+1}}.
\end{equation}
Eliminate $U$ using $q_r\, U_{a_r}^{[q_r]} = F_{a_r}^{[q_r-1]}$ to obtain an equivalent recursion written purely in terms of $F^{[\cdot]}$:
\begin{equation}\label{eq:F-recursion}
  \boxed{\;F_i^{[q]}
  = \sum_{m=1}^{q} \frac{1}{m!}
    \sum_{a_1=1}^d\!\cdots\!\sum_{a_m=1}^d
    \frac{\partial^{m} f_i(u_0)}{\partial u_{a_1}\cdots\partial u_{a_m}}
    \sum_{\substack{q_1+\cdots+q_m = q \\
                    q_r \ge 1\ \text{for } r=1,\dots,m}}
    \prod_{r=1}^{m} \frac{1}{q_r}\, F_{a_r}^{[q_r-1]}\;}\quad(q\ge 1),\qquad F_i^{[0]} = f_i(u_0).
\end{equation}
Remark (why the indices $a_r$ appear and why cross-component dependence is necessary):

- The chain rule differentiates $f_i(u)$ with respect to each input variable $u_{a_r}$ that actually changes along $u(t)$, hence the partial derivative indices $a_1,\dots,a_m$ in $\partial^{m} f_i/\partial u_{a_1}\cdots\partial u_{a_m}$.

- The product rule multiplies these partials by time-derivatives of the corresponding input components. Since $u'_{a}(t) = f_{a}(u(t))$, the Taylor coefficients of each input component $U_{a}^{[q]}$ are tied to the coefficients $F_{a}^{[q-1]}$ of the $a$-th right-hand side via \eqref{eq:U-F}. This is why each child factor is $\tfrac{1}{q_r} F_{a_r}^{[q_r-1]}$ and not something involving only index $i$.

- Consequently, $F_i^{[q]}$ depends on all components $a\in\{1,\dots,d\}$ through the mixed partials of $f_i$ and through the dynamics of each component $u_a$. It cannot depend only on $F_i^{[\cdot]}$ unless $f_i$ depends only on $u_i$.

\section{1.4 From the literal $F$-recursion to trees (complete constructive rewrite)}

We derive trees strictly from the literal recursion for $F^{[q]}$ and only then abstract step by step.

\section{1.4.1 Compute $F$ literally from Eq. \eqref{eq:F-recursion}}

Constraint handling: for each $m\in\{1,\dots,q\}$ enumerate all compositions $(q_1,\dots,q_m)$ of $q$ with $q_r\ge 1$ and $q_1+\cdots+q_m=q$.

Helper: `Compositions(q, m)` enumerates all $m$-tuples of positive integers summing to $q$.

```text
function Compute_F(k):  # returns table F[i][q] for i=1..d and q=0..k
  for i in {1,...,d}:            # base
    F[i][0] := f_i(u0)

  for q in {1,2,...,k}:          # increasing order
    for i in {1,...,d}:          # each output component
      S := 0
      for m in {1,...,q}:        # arity
        for (q1,...,qm) in Compositions(q, m):          # q1+...+qm=q, all q_r>=1
          for (a1,...,am) in {1,...,d}^m:               # index tuple
            layer := (1/m!) * (partial^m f_i(u0) / partial u_{a1}...partial u_{am})
            prod  := 1
            for r in {1,...,m}:
              prod := prod * (1/q_r) * F[a_r][q_r-1]
            S := S + layer * prod
      F[i][q] := S

  return F
```

This implements Eq. \eqref{eq:F-recursion} exactly, with the constraints $q_1+\cdots+q_m=q$ and $q_r\ge 1$ enforced by the `Compositions` loop.

\section{1.4.2 How rooted trees arise from the loop nest}

Direct expansion to $q=5$ (keep all $d$-sums explicit; lower orders unexpanded):
\begin{equation}
\label{eq:q5-direct}
\begin{aligned}
F_i^{[5]}
&= \underbrace{\frac{1}{5!} \sum_{a_1,\dots,a_5=1}^{d}
\frac{\partial^{5} f_i(u_0)}{\partial u_{a_1}\cdots\partial u_{a_5}}\,
F_{a_1}^{[0]} F_{a_2}^{[0]} F_{a_3}^{[0]} F_{a_4}^{[0]} F_{a_5}^{[0]}}_{m=5,\ (1,1,1,1,1)}
\\[4pt]
&\quad+ \underbrace{\frac{1}{4!} \sum_{a_1,\dots,a_4=1}^{d} \frac{\partial^{4} f_i(u_0)}{\partial u_{a_1}\cdots\partial u_{a_4}}\, \frac{1}{2}
\Bigl[ F_{a_1}^{[1]} F_{a_2}^{[0]} F_{a_3}^{[0]} F_{a_4}^{[0]}
+ F_{a_1}^{[0]} F_{a_2}^{[1]} F_{a_3}^{[0]} F_{a_4}^{[0]}
+ F_{a_1}^{[0]} F_{a_2}^{[0]} F_{a_3}^{[1]} F_{a_4}^{[0]}
+ F_{a_1}^{[0]} F_{a_2}^{[0]} F_{a_3}^{[0]} F_{a_4}^{[1]} \Bigr]}_{m=4,\ \text{perms of }(2,1,1,1)}
\\[4pt]
&\quad+ \underbrace{\frac{1}{3!} \sum_{a_1,a_2,a_3=1}^{d} \frac{\partial^{3} f_i(u_0)}{\partial u_{a_1}\partial u_{a_2}\partial u_{a_3}} \Bigl[
\tfrac{1}{3} \bigl( F_{a_1}^{[2]} F_{a_2}^{[0]} F_{a_3}^{[0]} + F_{a_1}^{[0]} F_{a_2}^{[2]} F_{a_3}^{[0]} + F_{a_1}^{[0]} F_{a_2}^{[0]} F_{a_3}^{[2]} \bigr)
+ \tfrac{1}{4} \bigl( F_{a_1}^{[1]} F_{a_2}^{[1]} F_{a_3}^{[0]} + F_{a_1}^{[1]} F_{a_2}^{[0]} F_{a_3}^{[1]} + F_{a_1}^{[0]} F_{a_2}^{[1]} F_{a_3}^{[1]} \bigr)
\Bigr]}_{m=3,\ (3,1,1)\ \text{and}\ (2,2,1)}
\\[4pt]
&\quad+ \underbrace{\frac{1}{2!} \sum_{a_1,a_2=1}^{d} \frac{\partial^{2} f_i(u_0)}{\partial u_{a_1}\partial u_{a_2}} \Bigl[
\tfrac{1}{4} \bigl( F_{a_1}^{[3]} F_{a_2}^{[0]} + F_{a_1}^{[0]} F_{a_2}^{[3]} \bigr)
+ \tfrac{1}{6} \bigl( F_{a_1}^{[2]} F_{a_2}^{[1]} + F_{a_1}^{[1]} F_{a_2}^{[2]} \bigr)
\Bigr]}_{m=2,\ (4,1)\ \text{and}\ (3,2)}
\\[4pt]
&\quad+ \underbrace{\sum_{a=1}^{d} \frac{\partial f_i(u_0)}{\partial u_{a}} \, \tfrac{1}{5} \, F_{a}^{[4]}}_{m=1,\ (5)}.
\end{aligned}
\end{equation}
Here $F^{[0]}=f(u_0)$; $F^{[1]}$, $F^{[2]}$, $F^{[3]}$, $F^{[4]}$ remain unexpanded. For each $m\in\{1,\dots,5\}$, we sum over all compositions $(q_1,\dots,q_m)$ of $5$ and over component indices; each child contributes $\tfrac{1}{q_r}\,F_{a_r}^{[q_r-1]}$ and the parent contributes $\tfrac{1}{m!}$.

Succinct interpretation.
- The five groups above correspond to $m=5,4,3,2,1$ with child-order patterns $(1,1,1,1,1)$, $(2,1,1,1)$, $(3,1,1)$ and $(2,2,1)$, $(4,1)$ and $(3,2)$, and $(5)$.
- This is the literal recurrence \eqref{eq:F-recursion}: higher-order coefficients are obtained by differentiating $f_i$ and attaching lower-order $F^{[\cdot]}$ with weights $\tfrac{1}{m!}\prod_r \tfrac{1}{q_r}$.

Rooted-tree shorthand (introduced after the expansion).
- Notation: write a node with $m$ children as $[\tau_1,\dots,\tau_m]$ and a leaf as $\bullet$ ($F^{[0]}$). Define $\operatorname{ord}(\bullet)=1$ and $\operatorname{ord}([\tau_1,\dots,\tau_m])=1+\sum_r \operatorname{ord}(\tau_r)$.
- Combinatorics (step-by-step):
  1) Fix $q\ge 1$ and consider trees $\tau$ with $\operatorname{ord}(\tau)=q+1$.
  2) Let the root have $m$ children $\tau_1,\dots,\tau_m$ with $\operatorname{ord}(\tau_r)\ge 1$. By definition,
     $\operatorname{ord}([\tau_1,\dots,\tau_m]) = 1 + \sum_{r=1}^m \operatorname{ord}(\tau_r)$, hence $\sum_{r=1}^m \operatorname{ord}(\tau_r) = q$.
  3) The ordered $m$-tuple $(\operatorname{ord}(\tau_1),\dots,\operatorname{ord}(\tau_m))$ is a composition of $q$ into $m$ positive integers. If we ignore the order of the children, this collapses to an integer partition of $q$.
  4) For $q=5$, the partitions are $1+1+1+1+1$, $2+1+1+1$, $3+1+1$, $2+2+1$, $4+1$, $3+2$, and $5$. Our loops enumerate all ordered versions (compositions) of these, and the single $1/m!$ factor at the root symmetrizes over permutations.

- Root child-order patterns at the root (unlabeled; children unordered):
  - $(1,1,1,1,1)$: $[\bullet,\bullet,\bullet,\bullet,\bullet]$.
  - $(2,1,1,1)$: $[[\bullet],\bullet,\bullet,\bullet]$.
  - $(3,1,1)$: a root with one child of order 3 and two leaves, e.g. $[[\bullet,\bullet],\bullet,\bullet]$.
  - $(2,2,1)$: two children of order 2 and one leaf, e.g. $[[\bullet],[\bullet],\bullet]$.
  - $(4,1)$: one child of order 4 and one leaf, e.g. $[[\bullet,\bullet,\bullet],\bullet]$.
  - $(3,2)$: one child of order 3 and one child of order 2, e.g. $[[\bullet,\bullet],[\bullet]]$.
  - $(5)$: a single child $[\tau]$ with $\operatorname{ord}(\tau)=5$. In the unexpanded formula, this group is exactly $\sum_{a} (\partial f_i/\partial u_a)(u_0)\, (1/5)\, F_a^{[4]}$, and all internal order-5 subtree shapes are contained within $F^{[4]}$.

- Evaluation (purely shorthand for the sums): for a given shape, $F_i([\tau_1,\dots,\tau_m])(u_0)$ denotes
  $$
  \sum_{a_1,\dots,a_m=1}^{d} \Bigl( \tfrac{1}{m!} \, \frac{\partial^{m} f_i(u_0)}{\partial u_{a_1}\cdots\partial u_{a_m}} \Bigr)
  \prod_{r=1}^{m} \Bigl( \tfrac{1}{q_r} \, F_{a_r}(\tau_r)(u_0) \Bigr),
  $$
  where $q_r=\operatorname{ord}(\tau_r)$ and $F_{a_r}(\bullet)=F_{a_r}^{[0]}$. For example, the entire $m=3$ contribution
  $$
  \frac{1}{3!} \sum_{a_1,a_2,a_3} \frac{\partial^3 f_i(u_0)}{\partial u_{a_1}\partial u_{a_2}\partial u_{a_3}} \Bigl[ \tfrac{1}{3} (F_{a_1}^{[2]}F_{a_2}^{[0]}F_{a_3}^{[0]} + \cdots) + \tfrac{1}{4} (F_{a_1}^{[1]}F_{a_2}^{[1]}F_{a_3}^{[0]} + \cdots) \Bigr]
  $$
  is precisely the sum over the two shapes $[\tau,\bullet,\bullet]$ with $\operatorname{ord}(\tau)=3$ and $[\sigma,\sigma',\bullet]$ with $\operatorname{ord}(\sigma)=\operatorname{ord}(\sigma')=2$, via the rule above. The shorthand does not add assumptions; it only indexes and names the already-present sums and weights in \eqref{eq:F-recursion}.

- Exact correspondence to \eqref{eq:F-recursion}: by construction,
  $$
  F_i^{[q]} \;=\; \sum_{m=1}^{q} \; \sum_{\substack{q_1+\cdots+q_m=q\\ q_r\ge 1}} F_i\bigl([\tau_1,\dots,\tau_m]\bigr)(u_0),\quad \text{with } \operatorname{ord}(\tau_r)=q_r.
  $$
  This is just a renaming of the loops: the outer $1/m!$ and inner $1/q_r$ factors are exactly those in \eqref{eq:F-recursion}.

Tree-form expansion of $F_i^{[5]}$ (shorthand that exactly equals the explicit sum in \eqref{eq:q5-direct}):
\begin{equation}
\label{eq:q5-tree}
\begin{aligned}
F_i^{[5]}
&= F_i([\bullet,\bullet,\bullet,\bullet,\bullet])(u_0)
\\[2pt]
&\quad+ F_i([[\bullet],\bullet,\bullet,\bullet])(u_0)
\\[2pt]
&\quad+ \Bigl( F_i([[\bullet,\bullet],\bullet,\bullet])(u_0) + F_i([[[\bullet]],\bullet,\bullet])(u_0) \Bigr)
\\[2pt]
&\quad+ F_i([[\bullet],[\bullet],\bullet])(u_0)
\\[2pt]
&\quad+ \Bigl( F_i([[\bullet,\bullet,\bullet],\bullet])(u_0) + F_i([[[\bullet],\bullet],\bullet])(u_0) + F_i([[[\bullet,\bullet]],\bullet])(u_0) + F_i([[[[\bullet]]],\bullet])(u_0) \Bigr)
\\[2pt]
&\quad+ \Bigl( F_i([[\bullet,\bullet],[\bullet]])(u_0) + F_i([[[\bullet]],[\bullet]])(u_0) \Bigr)
\\[2pt]
&\quad+ F_i([\tau])(u_0)\quad\text{with } \operatorname{ord}(\tau)=5\text{ (this equals the }m=1\text{ term } \sum_a (\partial f_i/\partial u_a)\, (1/5)\, F_a^{[4]}\text{)}.
\end{aligned}
\end{equation}

Python: rooted execution trees for the scheduling above
