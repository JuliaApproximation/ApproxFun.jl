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
$$u + {\rm e}^x \int_{-1}^1 \cos x u(x) {\rm d}x = \cos {\rm e}^x$$
We can solve this equation as follows:
```jldoctest
julia> Σ = DefiniteIntegral(Chebyshev()); x=Fun();

julia> u = (I+exp(x)*Σ[cos(x)])\cos(exp(x));

julia> u(0.1)
0.2186429485562879
```

## Boundary conditions

Incorporating boundary conditions into differential equations is important
so that the equation is well-posed.  This is accomplished via combining
operators and _functionals_ (i.e., `1 × ∞` operators).  As a simple example, consider the first order
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
6.455495187177958e-17
```
Behind the scenes, the `Vector{Operator{T}}` representing the Functionals
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
0.04999999999996024
```
Note in this case the space is inferred from the variable coefficient `x`.


## Systems of equations

Systems of equations can be handled by creating a matrix of operators and
functionals.  For example, we can solve the system

$$\begin{align*}
    u'' - u + 2v &= {\rm e}^x  \cr
    v' + v &= {\rm e}^x  \cr
    u(-1) &= u'(-1) = v(-1) = 0
\end{align*}$$
This can be solved via:
```jldoctest
julia> x = Fun(); B = Evaluation(Chebyshev(),-1);

julia> A = [B      0;
            B*𝒟    0;
            0      B;
            𝒟^2-I  2I;
            I      𝒟+I];

julia> u,v = A\[0;0;0;exp(x);cos(x)];

julia> u(-1),u'(-1),v(-1)
(6.938893903907228e-17,-2.7755575615628914e-16,0.0)

julia> norm(u''-u+2v-exp(x))
6.957393606436152e-16

julia> norm(u + v'+v-cos(x))
2.878538423331588e-16
```
In this example, the automatic space detection failed and so we needed
to specify the domain space for `B`.


## QR Factorization

Behind the scenes, `A\b` where `A` is an `Operator` is implemented via
an adaptive QR factorization.  That is, it is equivalent to
`qrfact(A)\b`.  Note that `qrfact` adaptively caches a partial QR Factorization
as it is applied to different right-hand sides, so the same operator can be
inverted much more efficiently in subsequent problems.


## Partial differential equations

Partial differential operators are also supported.  Here's an example
of solving the Poisson equation with zero boundary conditions:
```julia
d = (-1..1)^2
x,y = Fun(d)
f = exp.(-10(x+0.3)^2-20(y-0.2)^2)  # use broadcasting as exp(f) not implemented in 2D
A = [Dirichlet(d);Δ]              # Δ is an alias for Laplacian()
@time u = A \ [zeros(∂(d));f]     #4s for ~3k coefficients
```
This equation requires a tolerance to be solved in a reasonable amount of time
as there are weak singularities at the corners.  Using a QR Factorization
reduces the cost of subsequent calls substantially:
```julia
QR = qrfact(A)
@time QR \ [zeros(∂(d));f]   # 4s
g = exp.(-10(x+0.2)^2-20(y-0.1)^2)
@time \(QR,[zeros(∂(d));g];tolerance=1E-10)  # 0.09s
```
These timings are far from optimal, mostly due to the cost of Optimization of partial differential equations is currently work in progress.


## Nonlinear equations
