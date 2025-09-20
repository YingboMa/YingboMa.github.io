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

Python: rooted execution trees for the scheduling above
```python
from __future__ import annotations
from dataclasses import dataclass, field
from functools import lru_cache
from itertools import product
from typing import Tuple, List

# Unicode superscripts/subscripts for nicer derivative labels
_SUPERS = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")
_SUBS   = str.maketrans("0123456789-", "₀₁₂₃₄₅₆₇₈₉₋")

def sup_digits(n: int) -> str:
    return str(n).translate(_SUPERS)

def sub_digits(n: int) -> str:
    return str(n).translate(_SUBS)

@dataclass(frozen=True)
class TreeNode:
    kind: str  # 'leaf' | 'coeff' | 'deriv'
    label: str
    children: Tuple['TreeNode', ...] = field(default_factory=tuple)

    def pretty_lines(self, prefix: str = "", is_last: bool = True) -> List[str]:
        connector = "└─ " if is_last else "├─ "
        lines = [prefix + connector + self.label]
        new_prefix = prefix + ("   " if is_last else "│  ")
        for idx, child in enumerate(self.children):
            lines.extend(child.pretty_lines(new_prefix, idx == len(self.children) - 1))
        return lines

    def tau_summary(self) -> str:
        if self.kind == 'leaf':
            return "•"
        if self.kind == 'coeff':
            return self.children[0].tau_summary() if self.children else "•"
        # deriv node
        return "[" + ",".join(child.tau_summary() for child in self.children) + "]"

def compositions(n: int, m: int) -> List[Tuple[int, ...]]:
    if m == 0:
        return [()] if n == 0 else []
    result: List[Tuple[int, ...]] = []
    for first in range(1, n - m + 2):
        for rest in compositions(n - first, m - 1):
            result.append((first,) + rest)
    return result

@lru_cache(maxsize=None)
def enumerate_execution_trees(i: int, q: int, d: int) -> Tuple[TreeNode, ...]:
    # Base: F[i][0] = f_i(u0) corresponds to a single leaf
    if q == 0:
        return (TreeNode(kind='leaf', label=f"F[{i}][0] ≡ f{sub_digits(i)}(u{sub_digits(0)})"),)

    nodes: List[TreeNode] = []
    for m in range(1, q + 1):
        for q_tuple in compositions(q, m):
            for a_tuple in product(range(1, d + 1), repeat=m):
                child_options: List[Tuple[TreeNode, ...]] = []
                for r, q_r in enumerate(q_tuple):
                    a_r = a_tuple[r]
                    child_options.append(enumerate_execution_trees(a_r, q_r - 1, d))
                # Cartesian product over child expansions
                for chosen_children in product(*child_options):
                    # Wrap each child with its (1/q_r) factor to reflect the product rule weights
                    wrapped_children: Tuple[TreeNode, ...] = tuple(
                        TreeNode(kind='coeff', label=f"(1/{q_r}) ×", children=(chosen_children[idx],))
                        for idx, q_r in enumerate(q_tuple)
                    )
                    deriv = "∂" + sup_digits(m) + f" f{sub_digits(i)}/" + "".join(f"∂u{sub_digits(a)}" for a in a_tuple)
                    q_str = ",".join(str(x) for x in q_tuple)
                    node_label = f"(1/{m}!) {deriv}; q=({q_str})"
                    nodes.append(TreeNode(kind='deriv', label=node_label, children=wrapped_children))
    return tuple(nodes)

def print_execution_trees(i: int, q: int, d: int, limit: int | None = None) -> None:
    trees = list(enumerate_execution_trees(i, q, d))
    if limit is not None:
        trees = trees[:limit]
    print(f"Execution trees for F[{i}][{q}] with dimension d={d}: count={len(trees)}")
    for t in trees:
        print("τ =", t.tau_summary())
        for line in t.pretty_lines():
            print(line)

if __name__ == "__main__":
    # Small demonstration to keep output manageable
    print_execution_trees(i=1, q=2, d=2)
```
