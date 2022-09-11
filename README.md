# ApproxFun.jl

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaApproximation.github.io/ApproxFun.jl/stable)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://JuliaApproximation.github.io/ApproxFun.jl/latest)
[![Build Status](https://travis-ci.org/JuliaApproximation/ApproxFun.jl.svg?branch=master)](https://travis-ci.org/JuliaApproximation/ApproxFun.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/40qmoxp189pwtuda?svg=true)](https://ci.appveyor.com/project/dlfivefifty/approxfun-jl)
[![codecov](https://codecov.io/gh/JuliaApproximation/ApproxFun.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaApproximation/ApproxFun.jl)
[![deps](https://juliahub.com/docs/ApproxFun/deps.svg)](https://juliahub.com/ui/Packages/ApproxFun/jGqLz?t=2)
[![version](https://juliahub.com/docs/ApproxFun/version.svg)](https://juliahub.com/ui/Packages/ApproxFun/jGqLz)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![Join the chat at https://gitter.im/JuliaApproximation/ApproxFun.jl](https://badges.gitter.im/JuliaApproximation/ApproxFun.jl.svg)](https://gitter.im/JuliaApproximation/ApproxFun.jl?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)


ApproxFun is a package for approximating functions. It is in a similar vein to the Matlab
package [`Chebfun`](http://www.chebfun.org) and the Mathematica package [`RHPackage`](https://github.com/dlfivefifty/RHPackage).

The  [`ApproxFun Documentation`](https://JuliaApproximation.github.io/ApproxFun.jl/latest) contains detailed information, or read on for a brief overview of the package. The documentation contains examples of usage, such as solving ordinary and partial differential equations.

The  [`ApproxFun Examples`](https://github.com/JuliaApproximation/ApproxFunExamples) repo contains many examples of
using this package, in Jupyter notebooks and Julia scripts. Note that this is independently maintained, so it might not always be in sync with the latest version of `ApproxFun`. We recommend checking the examples in the documentation first, as these will always be compatible with the latest version of the package.

## Introduction

### Approximating Functions

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

plot(h; label="f + g^2")
scatter!(r, h.(r); label="roots")
scatter!(rp, h.(rp); label="extrema")
```

<img src=https://github.com/JuliaApproximation/ApproxFun.jl/raw/master/images/extrema.png width=500 height=400>


### Differentiation and integration


Notice from above that to find the extrema, we used `'` overridden for the `differentiate` function. Several other `Julia`
base functions are overridden for the purposes of calculus. We may check that the exponential is its own derivative, by evaluating the norm of the difference and checking that it is small:

```julia
f = Fun(exp, -1..1)
norm(f-f')  # 4.4391656415701095e-14
```

Similarly, `cumsum` defines an indefinite integration operator:

```julia
g = cumsum(f)
g = g + f(-1)
norm(f-g) # 3.4989733283850415e-15d
```

Algebraic and differential operations are also implemented where possible, and most of Julia's built-in functions (and special functions from [`SpecialFunctions.jl`](https://github.com/JuliaMath/SpecialFunctions.jl)) are overridden to accept `Fun`s:

```julia
x = Fun()
f = erf(x)
g = besselj(3,exp(f))
h = airyai(10asin(f)+2g)
```

## Examples of Usage

Check the [documentation](https://JuliaApproximation.github.io/ApproxFun.jl/latest) for examples of usage.

## References

J. L. Aurentz & R. M. Slevinsky (2019), On symmetrizing the ultraspherical spectral method for self-adjoint problems, arxiv:1903.08538

S. Olver & A. Townsend (2014), A practical framework for infinite-dimensional linear algebra, Proceedings of the 1st First Workshop for High Performance Technical Computing in Dynamic Languages, 57–62

A. Townsend & S. Olver (2014), The automatic solution of partial differential equations using a global spectral method,  J. Comp. Phys., 299: 106–123

S. Olver & A. Townsend (2013), Fast inverse transform sampling in one and two dimensions, arXiv:1307.1223

S. Olver & A. Townsend (2013), A fast and well-conditioned spectral method, SIAM Review, 55:462–489
