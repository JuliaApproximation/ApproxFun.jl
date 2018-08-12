# Linear equations

Linear equations such as ordinary and partial differential equations,
 fractional differential equations and integral equations can be solved using ApproxFun.
This is accomplished using `A\b` where `A` is an `Operator` and `b`
is a `Fun`.  As a simple example, consider the equation

$$u'(\theta) + cu(\theta) = \cos\theta$$

where we want a solution that is periodic on $[0,2\pi)$.  This can be solved succinctly
as follows:
```@meta
DocTestSetup = quote
    using ApproxFun
end
```
```jldoctest
julia> b = Fun(cos,Fourier());

julia> c = 0.1; u = (𝒟+c*I)\b;

julia> u(0.6)
0.64076835137228

julia> (c*cos(0.6)+sin(0.6))/(1+c^2)  # exact solution
0.6407683513722804
```
Recall that `𝒟` is an alias to `Derivative() == Derivative(UnsetSpace(),1)`.

As another example, consider the Fredholm integral equation

$$u + {\rm e}^x \int_{-1}^1 \cos x \, u(x) {\rm d}x = \cos {\rm e}^x$$

We can solve this equation as follows:
```jldoctest
julia> Σ = DefiniteIntegral(Chebyshev()); x=Fun();

julia> u = (I+exp(x)*Σ[cos(x)])\cos(exp(x));

julia> u(0.1)
0.21864294855628802
```
Note that we used the syntax `op[f::Fun]`, which is a shorthand for `op*Multiplication(f)`.

## Boundary conditions

Incorporating boundary conditions into differential equations is important
so that the equation is well-posed.  This is accomplished via combining
operators and functionals (i.e., `1 × ∞` operators).  As a simple example, consider the first order
initial value problem

$$u' = t u \qquad\hbox{and}\qquad u(0) = 1$$

To pose this in ApproxFun, we want to find a `u` such that `Evaluation(0)*u == 1`
and `(𝒟 - t)*u == 0`.  This is accomplished via:
```jldoctest
julia> t = Fun(0..1);

julia> u = [Evaluation(0); 𝒟 - t]  \ [1;0];

julia> u(0)
0.9999999999999996

julia> norm(u'-t*u)
1.2016080299388273e-16
```
Behind the scenes, the `Vector{Operator{T}}` representing the functionals
and operators are combined into a single `InterlaceOperator`.

A common usage is two-point boundary value problems. Consider
the singularly perturbed boundary value problem:

$$\epsilon u''-xu'+u = u \qquad u(-1) = 1,\quad u(1) = 2$$

This can be solved in ApproxFun via:
```jldoctest
julia> x = Fun();

julia> u = [Evaluation(-1);
            Evaluation(1);
            1/70*𝒟^2-x*𝒟+I] \ [1,2,0];

julia> u(0.1)
0.049999999999960326
```
Note in this case the space is inferred from the variable coefficient `x`.

This ODE can also be solved using the `Dirichlet` operator:
```jldoctest
julia> x = Fun();

julia> u = [Dirichlet();
            1/70*𝒟^2-x*𝒟+I] \ [[1,2],0];

julia> u(0.1)
0.04999999999996019
```

## Eigenvalue Problems

In analogy to linear algebra, many differential equations may be posed as eigenvalue problems. That is, for some differential operator $L$, there are a family of functions $u_i(x)$ such that
$$
L~u_i(x) = \lambda_i u_i(x)
$$
where $\lambda_i$ is the $i^{th}$ eigenvalue of the $L$ and has a corresponding *eigenfunction* $u_i(x)$. A classic eigenvalue problem is known as the quantum harmonic oscillator where
$$L = -\frac{1}{2}\frac{d^2}{dx^2} + \frac{1}{2} x^2$$
and one demands that $u(\infty) = u(-\infty) = 0$. Because we expect the solutions to be exponentially suppressed for large $x$, we can approximate this with Dirichlet boundary conditions at a 'reasonably large' $x$ without much difference.

We can express this in ApproxFun as the following:
```jldoctest
x = Fun(-8 .. 8)
L = -𝒟^2/2 + x^2/2
S = space(x)
B = Dirichlet(S)
λ, v = eigs(B, L, 500,tolerance=1E-10)
```
note that boundary conditions must be specified in the call to `eigs`. Plotting the first $20$ eigenfunctions offset vertically by their eigenvalue, we see

![harmonic_eigs](../assets/Harmonic_eigs.pdf)

If the solutions are not relatively constant near the boundary then one should push the boundaries further out.

For problems with different contraints or boundary conditions, `B` can be any zero functional constraint, eg. `DefiniteIntegral()`.

## Systems of equations

Systems of equations can be handled by creating a matrix of operators and
functionals.  For example, we can solve the system

$$\begin{align*}
    u'' - u + 2v &= {\rm e}^x  \cr
    v' + v &= cos(x) \cr
    u(-1) &= u'(-1) = v(-1) = 0
\end{align*}$$

using the following code:
```jldoctest
julia> x = Fun(); B = Evaluation(Chebyshev(),-1);

julia> A = [B      0;
            B*𝒟    0;
            0      B;
            𝒟^2-I  2I;
            I      𝒟+I];

julia> u,v = A\[0;0;0;exp(x);cos(x)];

julia> u(-1),u'(-1),v(-1)
(-4.163336342344337e-17,-2.7755575615628914e-16,-2.220446049250313e-16)

julia> norm(u''-u+2v-exp(x))
5.981056979045254e-16

julia> norm(u + v'+v-cos(x))
2.3189209621240424e-16
```
In this example, the automatic space detection failed and so we needed
to specify explicitly that the domain space for `B` is `Chebyshev()`.


## QR Factorization

Behind the scenes, `A\b` where `A` is an `Operator` is implemented via
an adaptive QR factorization.  That is, it is equivalent to
`qr(A)\b`.  (There is a subtly here in space inferring: `A\b` can use
    both `A` and `b` to determine the domain space, while `qr(A)` only
    sees the operator `A`.)
      Note that `qr` adaptively caches a partial QR Factorization
as it is applied to different right-hand sides, so the same operator can be
inverted much more efficiently in subsequent problems.


## Partial differential equations

Partial differential operators are also supported.  Here's an example
of solving the Poisson equation with zero boundary conditions:
```julia
d = Domain(-1..1)^2
x,y = Fun(d)
f = exp.(-10(x+0.3)^2-20(y-0.2)^2)  # use broadcasting as exp(f) not implemented in 2D
A = [Dirichlet(d);Δ]              # Δ is an alias for Laplacian()
@time u = A \ [zeros(∂(d));f]     #4s for ~3k coefficients
```
Using a QR Factorization
reduces the cost of subsequent calls substantially:
```julia
QR = qr(A)
@time QR \ [zeros(∂(d));f]   # 4s
g = exp.(-10(x+0.2)^2-20(y-0.1)^2)
@time QR \ [zeros(∂(d));g]  # 0.09s
```

Many PDEs have weak singularities at the corners, in which case it is beneficial to
specify a tolerance to reduce the time:
```julia
\(A,[zeros(∂(d));f]; tolerance=1E-6)
```


## Nonlinear equations

There is preliminary support for nonlinear equations, via Newton iteration
in function space.  Here is a simple two-point boundary value problem:

$$\begin{align*}
    \epsilon u'' &+ 6(1-x^2)u' +u^2=1 \cr
    u(-1)&=u(1)=0
\end{align*}$$

This can be solved using
```julia
x = Fun()
N = u->[u(-1.)-c;u(1.);ε*u''+6*(1-x^2)*u'+u^2-1.]
u = newton(N,u0)
```
