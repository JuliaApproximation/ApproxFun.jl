```@meta
DocTestSetup  = quote
    using ApproxFun, LinearAlgebra
end
```

# Operators

Linear operators between two spaces in ApproxFun are represented by subtypes of `Operator`.  Every operator has a `domainspace` and `rangespace`.  That is, if a `Fun` `f` has the space `domainspace(op)`, then`op*f` is a `Fun` with space `rangespace(op)`.

Note that the size of an operator is specified by the dimension of the domain and range space.

## Calculus operators

Differential and integral operators are perhaps the most useful type of operators in mathematics.  Consider the derivative operator on `CosSpace`:

```jldoctest def-D
julia> D = Derivative(CosSpace())
ConcreteDerivative : CosSpace(【0.0,6.283185307179586❫) → SinSpace(【0.0,6.283185307179586❫)
 ⋅  -1.0    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   ⋅
 ⋅    ⋅   -2.0    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   ⋅
 ⋅    ⋅     ⋅   -3.0    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   ⋅
 ⋅    ⋅     ⋅     ⋅   -4.0    ⋅     ⋅     ⋅     ⋅     ⋅   ⋅
 ⋅    ⋅     ⋅     ⋅     ⋅   -5.0    ⋅     ⋅     ⋅     ⋅   ⋅
 ⋅    ⋅     ⋅     ⋅     ⋅     ⋅   -6.0    ⋅     ⋅     ⋅   ⋅
 ⋅    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   -7.0    ⋅     ⋅   ⋅
 ⋅    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   -8.0    ⋅   ⋅
 ⋅    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   -9.0  ⋅
 ⋅    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   ⋱
 ⋅    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   ⋱

julia> f = Fun(θ->cos(cos(θ)), CosSpace());

julia> fp = D*f;

julia> fp(0.1) ≈ f'(0.1) ≈ sin(cos(0.1))*sin(0.1)
true
```

Here, we specified the domain space for the derivative operator, and it automatically
determined the range space:

```jldoctest def-D
julia> rangespace(D) == space(fp) == SinSpace()
true
```

Operators can be identified with infinite-dimensional matrices, whose entries are given by the canonical bases in the domain and range space.  In this case, the relevant formula is

```math
\mathop{D} \cos{kθ} = -k \sin{kθ}.
```

That is, the `(k,k+1)`th entry is as follows:

```jldoctest def-D
julia> k,j = 5,6;

julia> ej = Fun(domainspace(D), [zeros(j-1);1]);

julia> D[k,j] ≈ coefficient(D*ej, k) ≈ -k
true
```

The `Chebyshev` space has the property that its derivatives are given by ultraspherical spaces:

```jldoctest
julia> Derivative(Chebyshev())
ConcreteDerivative : Chebyshev() → Ultraspherical(1)
 ⋅  1.0   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   ⋅
 ⋅   ⋅   2.0   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   ⋅
 ⋅   ⋅    ⋅   3.0   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   ⋅
 ⋅   ⋅    ⋅    ⋅   4.0   ⋅    ⋅    ⋅    ⋅    ⋅   ⋅
 ⋅   ⋅    ⋅    ⋅    ⋅   5.0   ⋅    ⋅    ⋅    ⋅   ⋅
 ⋅   ⋅    ⋅    ⋅    ⋅    ⋅   6.0   ⋅    ⋅    ⋅   ⋅
 ⋅   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   7.0   ⋅    ⋅   ⋅
 ⋅   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   8.0   ⋅   ⋅
 ⋅   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   9.0  ⋅
 ⋅   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   ⋱
 ⋅   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   ⋱
```

## Functionals

A particularly useful class of operators are _functionals_, which map from functions to scalar numbers.  These are represented by operators of size `1 × ∞`: that is, infinite-dimensional analogues of row vectors.

As an example, the evaluation functional `f(0)` on `CosSpace` has the form:

```jldoctest def-D
julia> B = Evaluation(CosSpace(), 0)
ConcreteEvaluation : CosSpace(【0.0,6.283185307179586❫) → ConstantSpace(Point(0))
 1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  ⋯

julia> B*f ≈ f(0)
true
```

As can be seen from the output, `rangespace(B)` is a `ConstantSpace(Point(0))`, a one-dimensional space used to represent scalars whose domain is a single point, `0`.

Closely related to functionals are operators with finite-dimensional range.  For example, the `Dirichlet` operator represents the restriction of a space to its boundary.  In the case, of `Chebyshev()`, this amounts to evaluation at the endpoints `±1`:

```jldocetst
julia> B = Dirichlet(Chebyshev())
ConcreteDirichlet : Chebyshev() → 2-element ArraySpace:
ConstantSpace{Point{Float64}, Float64}[ConstantSpace(Point(-1.0)), ConstantSpace(Point(1.0))]
 1.0  -1.0  1.0  -1.0  1.0  -1.0  1.0  -1.0  1.0  -1.0  ⋯
 1.0   1.0  1.0   1.0  1.0   1.0  1.0   1.0  1.0   1.0  ⋯

julia> size(B)
(2, ℵ₀)

julia> B * Fun(exp) ≈ Fun([exp(-1), exp(1)])
true
```

## Multiplication

A `Multiplication` operator sends a `Fun` to a `Fun` in the corresponding space by multiplying a given function. The `Multiplication` operators are presented in matrix form in `ApproxFun`.

```jldoctest
julia> x = Fun();

julia> M = Multiplication(1 + 2x + x^2, Chebyshev())
ConcreteMultiplication : Chebyshev() → Chebyshev()
 1.5  1.0   0.25   ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅    ⋅
 2.0  1.75  1.0   0.25   ⋅     ⋅     ⋅     ⋅     ⋅     ⋅    ⋅
 0.5  1.0   1.5   1.0   0.25   ⋅     ⋅     ⋅     ⋅     ⋅    ⋅
  ⋅   0.25  1.0   1.5   1.0   0.25   ⋅     ⋅     ⋅     ⋅    ⋅
  ⋅    ⋅    0.25  1.0   1.5   1.0   0.25   ⋅     ⋅     ⋅    ⋅
  ⋅    ⋅     ⋅    0.25  1.0   1.5   1.0   0.25   ⋅     ⋅    ⋅
  ⋅    ⋅     ⋅     ⋅    0.25  1.0   1.5   1.0   0.25   ⋅    ⋅
  ⋅    ⋅     ⋅     ⋅     ⋅    0.25  1.0   1.5   1.0   0.25  ⋅
  ⋅    ⋅     ⋅     ⋅     ⋅     ⋅    0.25  1.0   1.5   1.0   ⋱
  ⋅    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅    0.25  1.0   1.5   ⋱
  ⋅    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋱     ⋱    ⋱

julia> coefficients(M * x) == coefficients((1 + 2x + x^2) * x) == M[1:4,1:2] * coefficients(x)
true
```

It is possible for domain space and range space to be different under `Mulitplication`.

```jldoctest
julia> cosθ = Fun(cos, CosSpace());

julia> sinθ = Fun(sin, SinSpace());

julia> sin2θ = Fun(x->sin(2x), SinSpace());

julia> cosθM = Multiplication(cosθ, SinSpace());

julia> cosθM * 2sinθ ≈ sin2θ
true
```

If a function is given by the expansion

```math
\mathop{f}(θ) = \sum_{n=1}^∞  f_n \sin{nθ}.
```

Then the matrix above can be easily derived from

```math
\begin{aligned}
\cos{θ} \mathop{f}(θ) &= \cos{θ} \sum_{n=1}^∞ f_n \sin{nθ} \\
                          &= \sum_{n=1}^∞ f_n \cos{θ} \sin{nθ} \\
                          &= \sum_{n=1}^∞ \frac{1}{2} f_n \left(\sin{(n-1)θ} + \sin{(n+1)θ}\right) \\
                          &= \sum_{n=1}^∞ \frac{1}{2} \left(f_{n-1} + f_{n+1}\right) \sin{nθ},
\end{aligned}
```

where ``f_0 = 0``.

## Algebraic manipulation of operators

Operators can be algebraically manipulated, provided that the domain and range spaces are compatible, or can be made compatible.  As a simple example, we can add the second derivative of a Fourier space to the
identity operator:

```jldoctest
julia> D2 = Derivative(Fourier(),2)
DerivativeWrapper : Fourier(【0.0,6.283185307179586❫) → Fourier(【0.0,6.283185307179586❫)
 0.0    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅      ⋅      ⋅      ⋅   ⋅
  ⋅   -1.0    ⋅     ⋅     ⋅     ⋅     ⋅      ⋅      ⋅      ⋅   ⋅
  ⋅     ⋅   -1.0    ⋅     ⋅     ⋅     ⋅      ⋅      ⋅      ⋅   ⋅
  ⋅     ⋅     ⋅   -4.0    ⋅     ⋅     ⋅      ⋅      ⋅      ⋅   ⋅
  ⋅     ⋅     ⋅     ⋅   -4.0    ⋅     ⋅      ⋅      ⋅      ⋅   ⋅
  ⋅     ⋅     ⋅     ⋅     ⋅   -9.0    ⋅      ⋅      ⋅      ⋅   ⋅
  ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   -9.0     ⋅      ⋅      ⋅   ⋅
  ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   -16.0     ⋅      ⋅   ⋅
  ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅      ⋅   -16.0     ⋅   ⋅
  ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅      ⋅      ⋅   -25.0  ⋅
  ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅      ⋅      ⋅      ⋅   ⋱

julia> D2 + I
PlusOperator : Fourier(【0.0,6.283185307179586❫) → Fourier(【0.0,6.283185307179586❫)
 1.0   ⋅    ⋅     ⋅     ⋅     ⋅     ⋅      ⋅      ⋅      ⋅   ⋅
  ⋅   0.0   ⋅     ⋅     ⋅     ⋅     ⋅      ⋅      ⋅      ⋅   ⋅
  ⋅    ⋅   0.0    ⋅     ⋅     ⋅     ⋅      ⋅      ⋅      ⋅   ⋅
  ⋅    ⋅    ⋅   -3.0    ⋅     ⋅     ⋅      ⋅      ⋅      ⋅   ⋅
  ⋅    ⋅    ⋅     ⋅   -3.0    ⋅     ⋅      ⋅      ⋅      ⋅   ⋅
  ⋅    ⋅    ⋅     ⋅     ⋅   -8.0    ⋅      ⋅      ⋅      ⋅   ⋅
  ⋅    ⋅    ⋅     ⋅     ⋅     ⋅   -8.0     ⋅      ⋅      ⋅   ⋅
  ⋅    ⋅    ⋅     ⋅     ⋅     ⋅     ⋅   -15.0     ⋅      ⋅   ⋅
  ⋅    ⋅    ⋅     ⋅     ⋅     ⋅     ⋅      ⋅   -15.0     ⋅   ⋅
  ⋅    ⋅    ⋅     ⋅     ⋅     ⋅     ⋅      ⋅      ⋅   -24.0  ⋅
  ⋅    ⋅    ⋅     ⋅     ⋅     ⋅     ⋅      ⋅      ⋅      ⋅   ⋱
```

When the domain and range space are not the same, the identity operator becomes a conversion operator.  That is, to represent `D+I` acting on the Chebyshev space, we would do the following:

```jldoctest
julia> D = Derivative(Chebyshev())
ConcreteDerivative : Chebyshev() → Ultraspherical(1)
 ⋅  1.0   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   ⋅
 ⋅   ⋅   2.0   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   ⋅
 ⋅   ⋅    ⋅   3.0   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   ⋅
 ⋅   ⋅    ⋅    ⋅   4.0   ⋅    ⋅    ⋅    ⋅    ⋅   ⋅
 ⋅   ⋅    ⋅    ⋅    ⋅   5.0   ⋅    ⋅    ⋅    ⋅   ⋅
 ⋅   ⋅    ⋅    ⋅    ⋅    ⋅   6.0   ⋅    ⋅    ⋅   ⋅
 ⋅   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   7.0   ⋅    ⋅   ⋅
 ⋅   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   8.0   ⋅   ⋅
 ⋅   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   9.0  ⋅
 ⋅   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   ⋱
 ⋅   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   ⋱

julia> C = Conversion(Chebyshev(), Ultraspherical(1))
ConcreteConversion : Chebyshev() → Ultraspherical(1)
 1.0  0.0  -0.5    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   ⋅
  ⋅   0.5   0.0  -0.5    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   ⋅
  ⋅    ⋅    0.5   0.0  -0.5    ⋅     ⋅     ⋅     ⋅     ⋅   ⋅
  ⋅    ⋅     ⋅    0.5   0.0  -0.5    ⋅     ⋅     ⋅     ⋅   ⋅
  ⋅    ⋅     ⋅     ⋅    0.5   0.0  -0.5    ⋅     ⋅     ⋅   ⋅
  ⋅    ⋅     ⋅     ⋅     ⋅    0.5   0.0  -0.5    ⋅     ⋅   ⋅
  ⋅    ⋅     ⋅     ⋅     ⋅     ⋅    0.5   0.0  -0.5    ⋅   ⋅
  ⋅    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅    0.5   0.0  -0.5  ⋅
  ⋅    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅    0.5   0.0  ⋱
  ⋅    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅    0.5  ⋱
  ⋅    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   ⋱

julia> D + I
PlusOperator : Chebyshev() → Ultraspherical(1)
 1.0  1.0  -0.5    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   ⋅
  ⋅   0.5   2.0  -0.5    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   ⋅
  ⋅    ⋅    0.5   3.0  -0.5    ⋅     ⋅     ⋅     ⋅     ⋅   ⋅
  ⋅    ⋅     ⋅    0.5   4.0  -0.5    ⋅     ⋅     ⋅     ⋅   ⋅
  ⋅    ⋅     ⋅     ⋅    0.5   5.0  -0.5    ⋅     ⋅     ⋅   ⋅
  ⋅    ⋅     ⋅     ⋅     ⋅    0.5   6.0  -0.5    ⋅     ⋅   ⋅
  ⋅    ⋅     ⋅     ⋅     ⋅     ⋅    0.5   7.0  -0.5    ⋅   ⋅
  ⋅    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅    0.5   8.0  -0.5  ⋅
  ⋅    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅    0.5   9.0  ⋱
  ⋅    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅    0.5  ⋱
  ⋅    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   ⋱

julia> (D + I)[1:10, 1:10] == (D + C)[1:10, 1:10]
true
```

ApproxFun can automatically determine the spaces, so if one writes `D + I` it will translate it to `D + C`.

Now consider the Fredholm integral operator of the second kind:

```math
\mathop{L} \mathop{u} = \mathop{u} + \mathop{e}^x \int_{-1}^1 \mathop{u}(x) \mathop{dx}.
```

We can construct this using

```jldoctest
julia> x = Fun();

julia> Σ = DefiniteIntegral(Chebyshev());

julia> L = I + exp(x)*Σ;

julia> u = cos(10x^2);

julia> (L*u)(0.1) ≈ u(0.1) + exp(0.1)*sum(u)
true
```

Note that `DefiniteIntegral` is a functional, i.e., a 1 × ∞ operator.  when multiplied on the left by a function, it automatically constructs the operator ``\mathop{e}^x \int_{-1}^1 \mathop{f}(x) \mathop{dx}`` via

```jldoctest
julia> x = Fun();

julia> Σ = DefiniteIntegral(Chebyshev());

julia> M = Multiplication(exp(x),ConstantSpace())
ConcreteMultiplication : ConstantSpace → Chebyshev()
 1.26607
 1.13032
 0.271495
 0.0443368
 0.00547424
 0.000542926
 4.49773e-5
 3.19844e-6
 1.99212e-7
 1.10368e-8
  ⋮

julia> M*Σ
TimesOperator : Chebyshev() → Chebyshev()
 2.53213     0.0  -0.844044     0.0  …  0.0  -0.0401926    0.0  ⋯
 2.26064     0.0  -0.753545     0.0     0.0  -0.0358831    0.0  ⋱
 0.542991    0.0  -0.180997     0.0     0.0  -0.0086189    0.0  ⋱
 0.0886737   0.0  -0.0295579    0.0     0.0  -0.00140752   0.0  ⋱
 0.0109485   0.0  -0.00364949   0.0     0.0  -0.000173785  0.0  ⋱
 0.00108585  0.0  -0.000361951  0.0  …  0.0  -1.72358e-5   0.0  ⋱
 8.99546e-5  0.0  -2.99849e-5   0.0     0.0  -1.42785e-6   0.0  ⋱
 6.39687e-6  0.0  -2.13229e-6   0.0     0.0  -1.01538e-7   0.0  ⋱
 3.98425e-7  0.0  -1.32808e-7   0.0     0.0  -6.32421e-9   0.0  ⋱
 2.20735e-8  0.0  -7.35785e-9   0.0     0.0  -3.50374e-10  0.0  ⋱
  ⋮           ⋱     ⋱            ⋱   …   ⋱     ⋱            ⋱   ⋱
```

!!! note
    `Σ*exp(x)` applies the operator to a function.  To construct the operator that first multiplies by `exp(x)`, use `Σ[exp(x)]`.  This is equivalent to `Σ*Multiplication(exp(x),Chebyshev())`.

## Operators and space promotion

It is often more convenient to not specify a space explicitly, but rather infer it when the operator is used.  For example, we can construct `Derivative()`, which has the alias `𝒟`, and represents the first derivative on any space:

```jldoctest
julia> f = Fun(cos, Chebyshev(0..1));

julia> (𝒟 * f)(0.1) ≈ -sin(0.1)
true

julia> f = Fun(cos, Fourier());

julia> (𝒟*f)(0.1) ≈ -sin(0.1)
true
```

Behind the scenes, `Derivative()` is equivalent to `Derivative(UnsetSpace(),1)`.  When multiplying a function `f`, the domain space is promoted before multiplying, that is, `Derivative()*f` is equivalent to `Derivative(space(f))*f`.

This promotion of the domain space happens even when operators have spaces attached.  This facilitates the following construction:

```jldoctest
julia> D = Derivative(Chebyshev());

julia> D^2
ConcreteDerivative : Chebyshev() → Ultraspherical(2)
 ⋅  ⋅  4.0   ⋅    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   ⋅
 ⋅  ⋅   ⋅   6.0   ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   ⋅
 ⋅  ⋅   ⋅    ⋅   8.0    ⋅     ⋅     ⋅     ⋅     ⋅   ⋅
 ⋅  ⋅   ⋅    ⋅    ⋅   10.0    ⋅     ⋅     ⋅     ⋅   ⋅
 ⋅  ⋅   ⋅    ⋅    ⋅     ⋅   12.0    ⋅     ⋅     ⋅   ⋅
 ⋅  ⋅   ⋅    ⋅    ⋅     ⋅     ⋅   14.0    ⋅     ⋅   ⋅
 ⋅  ⋅   ⋅    ⋅    ⋅     ⋅     ⋅     ⋅   16.0    ⋅   ⋅
 ⋅  ⋅   ⋅    ⋅    ⋅     ⋅     ⋅     ⋅     ⋅   18.0  ⋅
 ⋅  ⋅   ⋅    ⋅    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   ⋱
 ⋅  ⋅   ⋅    ⋅    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   ⋱
 ⋅  ⋅   ⋅    ⋅    ⋅     ⋅     ⋅     ⋅     ⋅     ⋅   ⋱
```

Note that `rangespace(D) ≠ Chebyshev()`, hence the operators are not compatible.  Therefore, it has thrown away its domain space, and thus this is equivalent to `Derivative(rangespace(D))*D`.

## Concatenating operators

The concatenation functions `vcat`, `hcat` and `hvcat` are overriden for operators to represent the resulting combined operator, now with a `rangespace` or `domainspace` that is an `ArraySpace`.
