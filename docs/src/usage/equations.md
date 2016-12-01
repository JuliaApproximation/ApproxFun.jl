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




## Nonlinear equations
