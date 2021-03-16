+++
title = "Taylor Series, Fourier Series, and Cauchy's Integral Formula"
hascode = true
date = Date(2018, 12, 17)
rss = "Using integration to compute high order derivatives of analytical functions."
+++

@def tags = ["math", "complex analysis"]


# Taylor Series, Fourier Series, and Cauchy's Integral Formula
\blogdate{Dec 17, 2018}


> The shortest and the best way between two truths of the real domain often
> passes through the imaginary one.
>
> --- Paul Painlevé

A Taylor series is a linear combination of monomials
$$
f(z) = \sum_{i=1}^n c_{i} z^i.
$$

A Fourier series is a linear combination of sinusoidal functions
$$
g(\theta) = \sum_{i=1}^n c_{i} e^{i \theta n}.
$$

Note that if we restrict the complex variable $z$ to the unit circle
$e^{i\theta}$, we see a Fourier series in disguise. Namely, we have
$f(\theta) = g(e^{i\theta})$. Hence, real Taylor and Fourier series are
special cases of a complex Taylor series.

The coefficients of a Taylor series can be obtained by differentiating the
function $f$ at a point $a$:

$$
c_{n} = \frac{f^{(n)}(a)}{n!}.
$$

The coefficients of a Fourier series are obtained by integrating $g$ with a
sinusoidal wave that oscillates $n$ times:

$$
c_{n} = \frac{1}{2\pi}\int_{0}^{2\pi} g(\theta) e^{-i\theta n}\; d\theta.
$$

They must be equal, if we take $z$ to be on the unit circle. We now have:

$$
\frac{f^{(n)}(a)}{n!} = \frac{1}{2\pi}\int_{0}^{2\pi} g(\theta)
e^{-i\theta n}\; d\theta.
$$

Remarkably, this equation says that we can use integration to compute high
order derivatives of an analytical function.

We can rewrite the above equation's right-hand side in the form of a complex
path integral. To be precise, the contour is a circle which centered at $a \in
\mathbb{C}$ with radius $1$.

\begin{align}
z &= e^{i\theta} + a \\
dz &= ie^{i\theta}d\theta \\
\frac{1}{i(z-a)} &= d\theta \\
\frac{1}{2\pi}\int_{0}^{2\pi} g(\theta) e^{-i\theta n}\; d\theta &=
\frac{1}{2\pi i}\oint \frac{f(z)}{(z-a)^{n+1}} \; dz.
\end{align}

Interestingly, we have Cauchy's integral formula:

$$
f^{(n)}(a) = \frac{n!}{2\pi i} \oint \frac{f(z)}{(z-a)^{n+1}}\; dz.
$$

Oh, what a coincidence!?

To consolidate and conclude this remark, I want to demonstrate a numerical
example of it.

```julia
julia> using FFTW

julia> z = cis.(range(0, 2pi, length=2^8+1)[1:end-1]);

julia> real( fft(exp.(z))./length(z) )[1:7] # Taylor expansion via FFT
7-element Array{Float64,1}:
 1.0
 1.0
 0.5
 0.16666666666666674
 0.041666666666666685
 0.00833333333333336
 0.0013888888888889245

julia> map(n->inv(factorial(n)), 0:6)
7-element Array{Float64,1}:
 1.0
 1.0
 0.5
 0.16666666666666666
 0.041666666666666664
 0.008333333333333333
 0.001388888888888889
```
