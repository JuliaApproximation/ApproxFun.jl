```@meta
DocTestSetup  = quote
    using ApproxFun, LinearAlgebra
end
```

# Spaces

A `Space` is an abstract type whose subtypes indicate which space a function lives in. This typically corresponds to the span of a (possibly infinite) basis.

## Classical orthogonal polynomial spaces

`Chebyshev`, `Ultraspherical`, `Jacobi`, `Hermite`, and `Laguerre` are spaces corresponding to expansion in classical orthogonal polynomials.

Note that we always use the classical normalization: the basis are _not_ orthonormal.  This is because this normalization leads to rational recurrence relationships, which are more efficient than their normalized counterparts. See the [Digital Library of Mathematical Functions](https://dlmf.nist.gov/18) for more information.

### Chebyshev space

The default space in ApproxFun is `Chebyshev`, which represents expansions in Chebyshev polynomials:

```math
\mathop{f}(x) = \sum_{k=0}^∞ f_k \mathop{T}_k(x),
```

where ``\mathop{T}_k(x) = \cos(k\arccos{x})``, which are orthogonal polynomials with respect to the weight

```math
\frac{1}{\sqrt{1-x^2}} \quad\text{for}\quad -1 ≤ x ≤ 1.
```

Note that there is an intrinsic link between `Chebyshev` and `CosSpace`:

```math
\mathop{g}(θ) = \mathop{f}(\cos{θ}) = \sum_{k=0}^∞ f_k \cos{kθ}.
```

In other words:

```jldoctest
julia> f = Fun(exp, Chebyshev());

julia> g = Fun(CosSpace(), coefficients(f));  # specify the coefficients directly

julia> f(cos(0.1)) ≈ exp(cos(0.1))
true

julia> g(0.1) ≈ exp(cos(0.1))
true
```

### Ultraspherical spaces

A key tool for solving differential equations are the ultraspherical spaces, encoded as `Ultraspherical(λ)` for `λ ≠ 0`,
which can be defined by the span of derivatives of Chebyshev polynomials, or alternatively as polynomials orthogonal with respect to the weight ``(1-x^2)^{λ - \frac{1}{2}}`` for ``-1 ≤ x ≤ 1``.

Note that `Ultraspherical(1)` corresponds to the Chebyshev basis of the second kind: ``\mathop{U}_k(x) = \frac{\sin((k+1)\arccos{x})}{\sin(\arccos{x})}``.  The relationship with Chebyshev polynomials follows from trigonemetric identities: ``\mathop{T}_k'(x) = k \mathop{U}_{k-1}(x)``.

Converting between ultraspherical polynomials (with integer orders) is extremely efficient: it requires ``\mathop{O}(n)`` operations, where ``n`` is the number of coefficients.

### Jacobi spaces

`Jacobi(b,a)` represents the space spanned by the Jacobi polynomials, which are orthogonal polynomials with respect to the weight

```math
(1+x)^b(1-x)^a
```

using the standard normalization.

## Fourier and Laurent spaces

There are several different spaces to represent functions on periodic domains, which are typically a `PeriodicSegment`, `Circle` or `PeriodicLine`.

`CosSpace` represents expansion in cosine series:

```math
\mathop{f}(θ) = \sum_{k=0}^∞ f_k \cos{kθ}.
```

`SinSpace` represents expansion in sine series:

```math
\mathop{f}(θ) = \sum_{k=0}^∞ f_k \sin{(k+1)θ}.
```

`Fourier` represents functions that are sums of sines and cosines.  Note that if a function has the form

```math
\mathop{f}(θ) = f_0 + \sum_{k=1}^∞ \left(f_k^\mathrm{s} \sin{kθ} + f_k^\mathrm{c} \cos{kθ}\right),
```

then the coefficients of the resulting `Fun` are ordered as ``[f_0, f_1^\mathrm{s}, f_1^\mathrm{c}, …]``.  For example:

```jldoctest
julia> f = Fun(Fourier(), [1,2,3,4]);

julia> f(0.1) ≈ 1 + 2sin(0.1) + 3cos(0.1) + 4sin(2*0.1)
true
```

`Taylor` represents expansion with only non-negative complex exponential terms:

```math
\mathop{f}(θ) = \sum_{k=0}^∞ f_k \mathop{e}^{ikθ}.
```

`Hardy{false}` represents expansion with only negative complex exponential terms:

```math
\mathop{f}(θ) = \sum_{k=0}^∞ f_k \mathop{e}^{-i(k+1)θ}.
```

`Laurent` represents functions that are sums of complex exponentials.  Note that if a function has the form

```math
\mathop{f}(θ) = \sum_{k=-∞}^∞ f_k \mathop{e}^{ikθ},
```

then the coefficients of the resulting `Fun` are order as ``[f_0, f_{-1}, f_1, …]``.  For example:

```jldoctest
julia> f = Fun(Laurent(), [1,2,3,4]);

julia> f(0.1) ≈ 1 + 2exp(-im*0.1) + 3exp(im*0.1) + 4exp(-2im*0.1)
true
```

## Modifier spaces

Some spaces are built out of other spaces:

### `JacobiWeight`

`JacobiWeight(β,α,space)`  weights `space`, which is typically `Chebyshev()` or `Jacobi(b,a)`, by a Jacobi weight `(1+x)^α*(1-x)^β`: in other words, if the basis for `space` is ``\mathop{ψ}_k(x)`` and the domain is the unit interval `-1..1`, then the basis for `JacobiWeight(β,α,space)` is ``(1+x)^α(1-x)^β \mathop{ψ}_k(x)``. If the domain is
not the unit interval, then the basis is determined by mapping back to the unit interval: that is, if ``\mathop{M}(x)`` is the map dictated by `tocanonical(space, x)`, where the canonical domain is the unit interval, then the basis is ``(1+\mathop{M}(x))^α(1-\mathop{M}(x))^β \mathop{ψ}_k(x)``. For example, if the domain is another interval `a..b`, then

```math
\mathop{M}(x) = \frac{2x-b-a}{b-a},
```

and the basis becomes

```math
\left(\frac{2}{b-a}\right)^{α+β}  (x-a)^α (b-x)^β \mathop{ψ}_k(x).
```

### `SumSpace`

`SumSpace((space_1,space_2,…,space_n))` represents the direct sum of the spaces, where evaluation is defined by adding up each component. A simple example is the following, showing that the coefficients are stored by interlacing:

```jldoctest
julia> x = Fun(identity, -1..1);

julia> f = cos(x-0.1)*sqrt(1-x^2) + exp(x);

julia> space(f)  # isa SumSpace
(1-x^2)^0.5[Chebyshev(-1..1)] ⊕ Chebyshev(-1..1)

julia> a, b = components(f);

julia> a(0.2) ≈ cos(0.2-0.1)*sqrt(1-0.2^2)
true

julia> b(0.2) ≈ exp(0.2)
true

julia> f(0.2) ≈ a(0.2) + b(0.2)
true

julia> norm(coefficients(f)[1:2:end] - coefficients(a))
0.0

julia> norm(coefficients(f)[2:2:end] - coefficients(b))
0.0
```

More complicated examples may interlace the coefficients using a different strategy.  Note that it is difficult to represent the first component of function ``\mathop{f}`` by a Chebyshev series because the derivatives of ``\mathop{f}`` at its boundaries blow up, whereas the derivative of a polynomial is a polynomial.

Note that `Fourier` and `Laurent` are currently implemented as `SumSpace`, but this may change in the future.

### `ArraySpace`

`ArraySpace(::Array{<:Space})` represents the direct sum of the spaces, where evaluation is defined in an array-wise manner.  A simple example is the following:

```jldoctest
julia> x = Fun(identity, -1..1);

julia> f = [exp(x); sqrt(1-x^2)*cos(x-0.1)];

julia> space(f)
2-element ArraySpace:
Space{ClosedInterval{Int64}, Float64}[Chebyshev(-1..1), (1-x^2)^0.5[Chebyshev(-1..1)]]

julia> a, b = components(f);

julia> norm(f(0.5) - [a(0.5); b(0.5)])
0.0

julia> norm(coefficients(f)[1:2:end] - coefficients(a))
0.0

julia> norm(coefficients(f)[2:2:end] - coefficients(b))
0.0
```

More complicated examples may interlace the coefficients using a different strategy.

### `TensorSpace`

`TensorSpace((space_1,space_2))` represents the tensor product of two spaces.
See the documentation of [`TensorSpace`](@ref) for more details on how the coefficients are interlaced.
Note that more than two spaces is only partially supported.

## Unset space

`UnsetSpace` is a special space that is used as a stand in when a space has not yet been determined, particularly by operators.
