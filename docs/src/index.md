# ApproxFun.jl Documentation



ApproxFun is a package for approximating and manipulating functions,
and for solving differential and integral equations.  

## Introduction

A basic approach of computational mathematics that ApproxFun exploits is expansion
in a basis

$$f(x) \approx \sum_{k=1}^n c_k \psi_k(x)$$

Some traditional examples of bases $\psi_1(x),\psi_2(x),\ldots$ are
1. Taylor series: $1,z,z^2,\ldots$
2. Fourier series (for periodic functions on `0..2π`): $1,\sin x, \cos x, \sin 2 x, \ldots$
3. Chebyshev series (for non-periodic functions on `-1..1`): $1,x,\cos 2 \hbox{acos}\, x, \cos 3 \hbox{acos}\, x, \ldots$

In ApproxFun, functions are represented by a `Fun` with two components: `space`,
which dictates the basis and `coefficients` which is a finite vector of coefficients.  Note that each `Fun` can have a different length vector of
coefficients, allowing for approximation of many different functions to high
accuracy.  

The approximation by a `Fun` can be determined by a variety of methods:

(1) Explicitly specifying the coefficients:
```julia
julia> f = Fun(Taylor(),[1,2,3]) # Represents 1 + 2z + 3z^2
Fun(Taylor(🕒),[1.0,2.0,3.0])

julia> f(1.0)
6.0
```
(2) Constructors take in a `Function` and adaptively determine the
    number of coefficients.  For example,
```julia
julia> Fun(exp)
Fun(Chebyshev(【-1.0,1.0】),[1.26607,1.13032,0.271495,0.0443368,0.00547424,0.000542926,4.49773e-5,3.19844e-6,1.99212e-7,1.10368e-8,5.5059e-10,2.49797e-11,1.03911e-12,3.99195e-14])
```
determines that `f` can be approximated to roughly machine precision using
14 coefficients.  See [Constructors](usage/constructors.md) for more information.

(3) Manipulation of `Fun`s give new `Fun`s, where the number of coefficients is determined from the input.  The simplest example is addition, which for compatible bases is just padding the vectors to the same length and adding.  
```julia
julia> a = Fun(cos,Chebyshev()); ncoefficients(a)
13

julia> b = Fun(x->cos(10cos(x^2)),Chebyshev()); ncoefficients(b)
51

julia> ncoefficients(a+b)
51
```
On the other hand, multiplication results in an approximation with more coefficients than either `a` or `b`, so that the result approximates the true `a*b` to roughly machine accuracy:
```
julia> ncoefficients(a*b)
63

julia> a(0.1)*b(0.1) - (a*b)(0.1)
1.1102230246251565e-16
```


The example of multiplication highlights the importance of adaptivity: if with a fixed discretization size, operations like multiplication would lose accuracy when the true function is no longer resolved by the discretization.  More complicated examples are solving differential equations, where the
coefficients of the solution can be determined adaptively, see [Equations](usage/equations.md).



ApproxFun supports a number of different spaces, as described in [Spaces](usage/spaces.md).  A key component of ApproxFun is support for interaction between different
spaces.  This is crucial for efficient solution of differential equations, where linear operators are described as acting between different spaces, see [Operators](usage/operators.md).  


## Contents



```@contents
Pages = ["usage/constructors.md",
         "usage/domains.md",
         "usage/spaces.md",
         "usage/operators.md",
         "usage/equations.md",
         "faq.md",
         "library.md"]
```
