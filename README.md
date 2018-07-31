# ApproxFun.jl

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaApproximation.github.io/ApproxFun.jl/stable)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://JuliaApproximation.github.io/ApproxFun.jl/latest)
[![Build Status](https://travis-ci.org/JuliaApproximation/ApproxFun.jl.svg?branch=master)](https://travis-ci.org/JuliaApproximation/ApproxFun.jl) [![Coverage Status](https://img.shields.io/coveralls/JuliaApproximation/ApproxFun.jl.svg)](https://coveralls.io/r/JuliaApproximation/ApproxFun.jl?branch=master) [![Join the chat at https://gitter.im/JuliaApproximation/ApproxFun.jl](https://badges.gitter.im/JuliaApproximation/ApproxFun.jl.svg)](https://gitter.im/JuliaApproximation/ApproxFun.jl?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)



ApproxFun is a package for approximating functions. It is in a similar vein to the Matlab
package [`Chebfun`](http://www.chebfun.org) and the Mathematica package [`RHPackage`](https://github.com/dlfivefifty/RHPackage).  

The  [`ApproxFun Documentation`](https://JuliaApproximation.github.io/ApproxFun.jl/latest) contains detailed information, or read on for a brief overview of the package.

The  [`ApproxFun Examples`](https://github.com/JuliaApproximation/ApproxFunExamples) contains many examples of
using this package, in Jupyter notebooks and Julia scripts.

## Introduction


Take your two favourite functions on an interval and create approximations to them as simply as:

```julia
using LinearAlgebra, SpecialFunctions, Plots, ApproxFun
x = Fun(identity,0..10)
f = sin(x^2)
g = cos(x)
```

Evaluating `f(.1)` will return a high
accuracy approximation to `sin(0.01)`. All the algebraic manipulations of functions
are supported and more.  For example, we can add `f` and `g^2` together and compute
the roots and extrema:

```julia
h = f + g^2
r = roots(h)
rp = roots(h')

using Plots
plot(h)
scatter!(r,h.(r))
scatter!(rp,h.(rp))
```

<img src=https://github.com/JuliaApproximation/ApproxFun.jl/raw/master/images/extrema.png width=500 height=400>


## Differentiation and integration


Notice from above that to find the extrema, we used `'` overridden for the `differentiate` function. Several other `Julia`
base functions are overridden for the purposes of calculus. Because the exponential is its own
derivative, the `norm` is small:

```julia
f = Fun(x->exp(x),-1..1)
norm(f-f')
```

Similarly, `cumsum` defines an indefinite integration operator:

```julia
g = cumsum(f)
g = g + f(-1)
norm(f-g)
```

Algebraic and differential operations are also implemented where possible, and most of Julia's built-in functions are overridden to accept `Fun`s:

```julia
x = Fun()
f = erf(x)
g = besselj(3,exp(f))
h = airyai(10asin(f)+2g)
```


## Solving ordinary differential equations


Solve the Airy ODE `u'' - x u = 0` as a BVP on `[-1000,200]`:

```julia
x = Fun(identity,-1000..200)
d = domain(x)
D = Derivative(d)
B = Dirichlet(d)
L = D^2 - x
u = [B;L] \ [airy.(endpoints(d)), 0]
plot(u)
```

<img src=https://github.com/JuliaApproximation/ApproxFun.jl/raw/master/images/airy.png width=500 height=400>


## Nonlinear Boundary Value problems

Solve a nonlinear boundary value problem satisfying the ODE `0.001u'' + 6*(1-x^2)*u' + u^2 = 1` with boundary conditions `u(-1)==1`, `u(1)==-0.5` on `[-1,1]`:

```julia
x=Fun()
u0=0.0x

N=u->[u(-1)-1,u(1)+0.5,0.001u''+6*(1-x^2)*u'+u^2-1]
u=newton(N,u0)
plot(u)
```

<img src=https://github.com/JuliaApproximation/ApproxFun.jl/raw/master/images/nbvp.png width=500 height=400>



## Periodic functions


There is also support for Fourier representations of functions on periodic intervals.
Specify the space `Fourier` to ensure that the representation is periodic:

```julia
f = Fun(cos,Fourier(-π..π))
norm(f' + Fun(sin,Fourier(-π..π))
```

Due to the periodicity, Fourier representations allow for the asymptotic savings of `2/π`
in the number of coefficients that need to be stored compared with a Chebyshev representation.
ODEs can also be solved when the solution is periodic:

```julia
s = Chebyshev(-π..π)
a = Fun(t-> 1+sin(cos(2t)),s)
L = Derivative() + a
f = Fun(t->exp(sin(10t)),s)
B = periodic(s,0)
uChebyshev = [B;L]\[0.;f]

s = Fourier(-π..π)
a = Fun(t-> 1+sin(cos(2t)),s)
L = Derivative() + a
f = Fun(t->exp(sin(10t)),s)
uFourier = L\f

ncoefficients(uFourier)/ncoefficients(uChebyshev),2/π
plot(uFourier)
```

<img src=https://github.com/JuliaApproximation/ApproxFun.jl/raw/master/images/periodic.png width=500 height=400>



## Sampling


Other operations including random number sampling using [Olver & Townsend 2013].  The
following code samples 10,000 from a PDF given as the absolute value of the sine function on `[-5,5]`:

```julia
f = abs(Fun(sin,-5..5))
x = ApproxFun.sample(f,10000)
histogram(x;normed=true)
plot!(f/sum(f))
```

<img src=https://github.com/JuliaApproximation/ApproxFun.jl/raw/master/images/sample.png width=500 height=400>


## Solving partial differential equations


We can solve PDEs, the following solves Helmholtz `Δu + 100u=0` with `u(±1,y)=u(x,±1)=1`
on a square.  This function has weak singularities at the corner,
so we specify a lower tolerance to avoid resolving these singularities
completely.

```julia
d = ChebyshevInterval()^2                            # Defines a rectangle
Δ = Laplacian(d)                            # Represent the Laplacian
f = ones(∂(d))                              # one at the boundary
u = \([Dirichlet(d);Δ+100I],[f;0.];         # Solve the PDE
                tolerance=1E-5)             
surface(u)                                  # Surface plot
```

<img src=https://github.com/JuliaApproximation/ApproxFun.jl/raw/master/images/helmholtz.png width=500 height=400>


<!-- We can also evolve PDEs.  The following solves advection—diffusion
`u_t = 0.01Δu - 4u_x -3u_y` on a rectangle

```julia
d = ChebyshevInterval()^2
u0 = Fun((x,y)->exp(-40(x-.1)^2-40(y+.2)^2),d)
B = dirichlet(d)
D = Derivative(ChebyshevInterval())
L = (0.01D^2-4D)⊗I + I⊗(0.01D^2-3D)
h = 0.002
timeevolution(B,L,u0,h)                    # Requires GLPlot
``` -->


## High precision

Solving differential equations with high precision types is available.  The following calculates `e` to 300 digits by solving the ODE `u' = u`:

```julia
setprecision(1000) do
    d=BigFloat(0)..BigFloat(1)
    D=Derivative(d)
    u=[ldirichlet();D-I]\[1;0]
    @test u(1) ≈ exp(BigFloat(1))
end
```



## References

S. Olver & A. Townsend (2014), A practical framework for infinite-dimensional linear algebra, Proceedings of the 1st First Workshop for High Performance Technical Computing in Dynamic Languages, 57–62

A. Townsend & S. Olver (2014), The automatic solution of partial differential equations using a global spectral method,  J. Comp. Phys., 299: 106–123

S. Olver & A. Townsend (2013), Fast inverse transform sampling in one and two dimensions, arXiv:1307.1223

S. Olver & A. Townsend (2013), A fast and well-conditioned spectral method, SIAM Review, 55:462–489
