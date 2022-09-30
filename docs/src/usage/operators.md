```@setup using-pkgs
using ApproxFun, LinearAlgebra
```

# Operators

Linear operators between two spaces in ApproxFun are represented by subtypes of `Operator`.  Every operator has a `domainspace` and `rangespace`.  That is, if a `Fun` `f` has the space `domainspace(op)`, then`op*f` is a `Fun` with space `rangespace(op)`.

Note that the size of an operator is specified by the dimension of the domain and range space.

## Calculus operators

Differential and integral operators are perhaps the most useful type of operators in mathematics.  Consider the derivative operator on `CosSpace`:

```@repl using-pkgs
D = Derivative(CosSpace())
f = Fun(θ->cos(cos(θ)), CosSpace());
fp = D*f;
fp(0.1) ≈ f'(0.1) ≈ sin(cos(0.1))*sin(0.1)
```

Here, we specified the domain space for the derivative operator, and it automatically
determined the range space:

```@setup def-D
using ApproxFun
D = Derivative(CosSpace())
f = Fun(θ->cos(cos(θ)),CosSpace())
fp = D*f
```

```@repl def-D
rangespace(D) == space(fp) == SinSpace()
```

Operators can be identified with infinite-dimensional matrices, whose entries are given by the canonical bases in the domain and range space.  In this case, the relevant formula is

```math
\mathop{D} \cos{kθ} = -k \sin{kθ}.
```

That is, the `(k,k+1)`th entry is as follows:

```@repl def-D
k,j = 5,6;
ej = Fun(domainspace(D),[zeros(j-1);1]);
D[k,j] ≈ (D*ej).coefficients[k] ≈ -k
```

The `Chebyshev` space has the property that its derivatives are given by ultraspherical spaces:

```@repl using-pkgs
Derivative(Chebyshev())
```

## Functionals

A particularly useful class of operators are _functionals_, which map from functions to scalar numbers.  These are represented by operators of size `1 × ∞`: that is, infinite-dimensional analogues of row vectors.

As an example, the evaluation functional `f(0)` on `CosSpace` has the form:

```@repl using-pkgs
B = Evaluation(CosSpace(),0)
B*f ≈ f(0)
```

As can be seen from the output, `rangespace(B)` is a `ConstantSpace(Point(0))`, a one-dimensional space used to represent scalars whose domain is a single point, `0`.

Closely related to functionals are operators with finite-dimensional range.  For example, the `Dirichlet` operator represents the restriction of a space to its boundary.  In the case, of `Chebyshev()`, this amounts to evaluation at the endpoints `±1`:

```@repl using-pkgs
B = Dirichlet(Chebyshev())
size(B)
B*Fun(exp)
B*Fun(exp) ≈ Fun([exp(-1),exp(1)])
```

## Multiplication

A `Multiplication` operator sends a `Fun` to a `Fun` in the corresponding space by multiplying a given function. The `Multiplication` operators are presented in matrix form in `ApproxFun`.

```@repl using-pkgs
x = Fun();
M = Multiplication(1 + 2x + x^2, Chebyshev())
(M * x).coefficients == ((1 + 2x + x^2) * x).coefficients == M[1:4,1:2] * x.coefficients
```

It is possible for domain space and range space to be different under `Mulitplication`.

```@repl using-pkgs
c = Fun(θ -> cos(θ), CosSpace());
Multiplication(c, SinSpace())
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

```@repl using-pkgs
D2 = Derivative(Fourier(),2)
D2 + I
```

When the domain and range space are not the same, the identity operator becomes a conversion operator.  That is, to represent `D+I` acting on the Chebyshev space, we would do the following:

```@repl using-pkgs
D = Derivative(Chebyshev())
C = Conversion(Chebyshev(),Ultraspherical(1))
D + C
```

ApproxFun can automatically determine the spaces, so if one writes `D + I` it will translate it to `D + C`.

Now consider the Fredholm integral operator of the second kind:

```math
\mathop{L} \mathop{u} = \mathop{u} + \mathop{e}^x \int_{-1}^1 \mathop{u}(x) \mathop{dx}.
```

We can construct this using

```@repl using-pkgs
x = Fun();
Σ = DefiniteIntegral(Chebyshev())
L = I + exp(x)*Σ
u = cos(10x^2);
(L*u)(0.1)
u(0.1) + exp(0.1)*sum(u)
```

Note that `DefiniteIntegral` is a functional, i.e., a 1 × ∞ operator.  when multiplied on the left by a function, it automatically constructs the operator ``\mathop{e}^x \int_{-1}^1 \mathop{f}(x) \mathop{dx}`` via

```@setup def-Σ
using ApproxFun
x = Fun()
Σ = DefiniteIntegral(Chebyshev())
```

```@repl def-Σ
M = Multiplication(exp(x),ConstantSpace())
M*Σ
```

Note that `Σ*exp(x)` applies the operator to a function.  To construct the operator that first multiplies by `exp(x)`, use `Σ[exp(x)]`.  This is equivalent to `Σ*Multiplication(exp(x),Chebyshev())`.

## Operators and space promotion

It is often more convenient to not specify a space explicitly, but rather infer it when the operator is used.  For example, we can construct `Derivative()`, which has the alias `𝒟`, and represents the first derivative on any space:

```@repl using-pkgs
f = Fun(cos,Chebyshev(0..1)); (𝒟*f)(0.1)
f = Fun(cos,Fourier()); (𝒟*f)(0.1)
```

Behind the scenes, `Derivative()` is equivalent to `Derivative(UnsetSpace(),1)`.  When multiplying a function `f`, the domain space is promoted before multiplying, that is, `Derivative()*f` is equivalent to `Derivative(space(f))*f`.

This promotion of the domain space happens even when operators have spaces attached.  This facilitates the following construction:

```jldoctest; setup=:(using ApproxFun)
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
