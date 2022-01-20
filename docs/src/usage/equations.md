```@setup using-pkgs
using ApproxFun, LinearAlgebra
```

# Linear equations

Linear equations such as ordinary and partial differential equations, fractional differential equations and integral equations can be solved using ApproxFun.  This is accomplished using `A\b` where `A` is an `Operator` and `b` is a `Fun`.  As a simple example, consider the equation

```math
\mathop{u}'(θ) + c \mathop{u}(θ) = \cos{θ},
```

where we want a solution that is periodic on ``[0,2π)``.  This can be solved succinctly as follows:

```@repl using-pkgs
b = Fun(cos,Fourier());
c = 0.1; u = (𝒟+c*I)\b;
u(0.6)
(c*cos(0.6)+sin(0.6))/(1+c^2)  # exact solution
```

Recall that `𝒟` is an alias to `Derivative() == Derivative(UnsetSpace(),1)`.

As another example, consider the Fredholm integral equation

```math
\mathop{u} + \mathop{e}^x \int_{-1}^1 \cos{x} \mathop{u}(x) \mathop{dx} = \cos{\mathop{e}^x}.
```

We can solve this equation as follows:

```@repl using-pkgs
Σ = DefiniteIntegral(Chebyshev()); x = Fun();
u = (I+exp(x)*Σ[cos(x)]) \ cos(exp(x));
u(0.1)
```

Note that we used the syntax `op[f::Fun]`, which is a shorthand for `op*Multiplication(f)`.

## Boundary conditions

Incorporating boundary conditions into differential equations is important so that the equation is well-posed.  This is accomplished via combining operators and functionals (i.e., `1 × ∞` operators).  As a simple example, consider the first order initial value problem

```math
\begin{gathered}
    \mathop{u}' = t \mathop{u}, \\
    \mathop{u}(0) = 1.
\end{gathered}
```

To pose this in ApproxFun, we want to find a `u` such that `Evaluation(0)*u == 1` and `(𝒟 - t)*u == 0`.  This is accomplished via:

```@repl using-pkgs
t = Fun(0..1);
u = [Evaluation(0); 𝒟 - t]  \ [1;0];
u(0)
norm(u'-t*u)
```

Behind the scenes, the `Vector{Operator{T}}` representing the functionals and operators are combined into a single `InterlaceOperator`.

A common usage is two-point boundary value problems. Consider the singularly perturbed boundary value problem:

```math
\begin{gathered}
    ϵ \mathop{u}'' - x \mathop{u}' + \mathop{u} = 0, \\
    \mathop{u}(-1) = 1, \mathop{u}(1) = 2.
\end{gathered}
```

This can be solved in ApproxFun via:

```@repl using-pkgs
x = Fun();
u = [Evaluation(-1);
     Evaluation(1);
     1/70*𝒟^2-x*𝒟+I] \ [1,2,0];
u(0.1)
```

Note in this case the space is inferred from the variable coefficient `x`.

This ODE can also be solved using the `Dirichlet` operator:

```@repl using-pkgs
x = Fun();
u = [Dirichlet();
     1/70*𝒟^2-x*𝒟+I] \ [[1,2],0];
u(0.1)
```

## Eigenvalue Problems

In analogy to linear algebra, many differential equations may be posed as eigenvalue problems. That is, for some differential operator ``\mathop{L}``, there are a family of functions ``\mathop{u}_i(x)`` such that

```math
\mathop{L} \mathop{u}_i(x) = λ_i \mathop{u}_i(x),
```

where ``λ_i`` is the ``i^{th}`` eigenvalue of the ``L`` and has a corresponding *eigenfunction* ``\mathop{u}_i(x)``. A classic eigenvalue problem is known as the quantum harmonic oscillator where

```math
\mathop{L} = -\frac{1}{2}\frac{\mathop{d}^2}{\mathop{dx}^2} + \frac{1}{2} x^2,
```

and one demands that ``\mathop{u}(∞) = \mathop{u}(-∞) = 0``. Because we expect the solutions to be exponentially suppressed for large ``x``, we can approximate this with Dirichlet boundary conditions at a 'reasonably large' ``x`` without much difference.

We can express this in ApproxFun as the following:

```julia
x = Fun(-8 .. 8)
L = -𝒟^2/2 + x^2/2
S = space(x)
B = Dirichlet(S)
λ, v = ApproxFun.eigs(B, L, 500,tolerance=1E-10)
```

Note that boundary conditions must be specified in the call to `eigs`. Plotting the first ``20`` eigenfunctions offset vertically by their eigenvalue, we see

![harmonic_eigs](../assets/Harmonic_eigs.pdf)

If the solutions are not relatively constant near the boundary then one should push the boundaries further out.

For problems with different contraints or boundary conditions, `B` can be any zero functional constraint, e.g., `DefiniteIntegral()`.

## Systems of equations

Systems of equations can be handled by creating a matrix of operators and functionals.  For example, we can solve the system

```math
\begin{gathered}
    \mathop{u}'' - \mathop{u} + 2 \mathop{v} = \mathop{e}^x,  \\
    \mathop{v}' + \mathop{v} = \cos{x}, \\
    \mathop{u}(-1) = \mathop{u}'(-1) = \mathop{v}(-1) = 0
\end{gathered}
```

using the following code:

```@repl using-pkgs
x = Fun(); B = Evaluation(Chebyshev(),-1);
A = [B      0;
     B*𝒟    0;
     0      B;
     𝒟^2-I  2I;
     I      𝒟+I];
u,v = A\[0;0;0;exp(x);cos(x)];
u(-1),u'(-1),v(-1)
norm(u''-u+2v-exp(x))
norm(u + v'+v-cos(x))
```

In this example, the automatic space detection failed and so we needed to specify explicitly that the domain space for `B` is `Chebyshev()`.

## QR Factorization

Behind the scenes, `A\b` where `A` is an `Operator` is implemented via an adaptive QR factorization.  That is, it is equivalent to `qr(A)\b`.  (There is a subtlety here in space inferring: `A\b` can use both `A` and `b` to determine the domain space, while `qr(A)` only sees the operator `A`.)  Note that `qr` adaptively caches a partial QR Factorization
as it is applied to different right-hand sides, so the same operator can be inverted much more efficiently in subsequent problems.

## Partial differential equations

Partial differential operators are also supported.  Here's an example
of solving the Poisson equation with zero boundary conditions:

```julia
d = Domain(-1..1)^2
x,y = Fun(d)
f = exp.(-10(x+0.3)^2-20(y-0.2)^2)  # use broadcasting as exp(f) not implemented in 2D
A = [Dirichlet(d);Δ]  # Δ is an alias for Laplacian()
@time u = A \ [zeros(∂(d));f]  # 4s for ~3k coefficients
```

Using a QR Factorization reduces the cost of subsequent calls substantially:

```julia
QR = qr(A)
@time QR \ [zeros(∂(d));f]  # 4s
g = exp.(-10(x+0.2)^2-20(y-0.1)^2)
@time QR \ [zeros(∂(d));g]  # 0.09s
```

Many PDEs have weak singularities at the corners, in which case it is beneficial to specify a tolerance to reduce the time:

```julia
\(A, [zeros(∂(d));f]; tolerance=1E-6)
```

## Nonlinear equations

There is preliminary support for nonlinear equations, via Newton iteration in function space.  Here is a simple two-point boundary value problem:

```math
\begin{gathered}
    ϵ \mathop{u}'' + 6(1-x^2) \mathop{u}' + \mathop{u}^2=1, \\
    \mathop{u}(-1) = \mathop{u}(1) = 0.
\end{gathered}
```

This can be solved using

```julia
x = Fun()
N = u -> [u(-1.)-c; u(1.); ε*u'' + 6*(1-x^2)*u' + u^2 - 1.0]
u = newton(N,u0)
```
