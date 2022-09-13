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
c = 0.1; u = (𝒟+c*I) \ b;
u(0.6)
(c*cos(0.6)+sin(0.6)) / (1+c^2)  # exact solution
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
u = [Evaluation(0); 𝒟-t] \ [1;0];
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
ϵ = 1/70; x = Fun();
u = [Evaluation(-1);
     Evaluation(1);
     ϵ*𝒟^2-x*𝒟+I] \ [1,2,0];
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

## QR Factorization

Behind the scenes, `A\b` where `A` is an `Operator` is implemented via an adaptive QR factorization.  That is, it is equivalent to `qr(A)\b`.  (There is a subtlety here in space inferring: `A\b` can use both `A` and `b` to determine the domain space, while `qr(A)` only sees the operator `A`.)  Note that `qr` adaptively caches a partial QR Factorization
as it is applied to different right-hand sides, so the same operator can be inverted much more efficiently in subsequent problems.
