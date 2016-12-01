# Spaces

A `Space` is an abstract type whose subtypes indicate which space a function lives in.
This typically corresponds to the span of a (possibly infinite) basis.

## Chebyshev space

The default space in ApproxFun is `Chebyshev`, which represents expansions
in Chebyshev polynomials:

$$f(x) = \sum_{k=0}^\infty f_k T_k(x)$$

where $T_k(x) = \cos k \,{\rm acos} x$.  
Note that there is an intrinsic link between `Chebyshev` and `CosSpace`:  

$$g(\theta) = f(\cos \theta) = \sum_{k=0}^\infty f_k \cos k \theta$$

In other words:
```@meta
DocTestSetup = quote
    using ApproxFun
end
```
```jldoctest
julia> f=Fun(exp,Chebyshev());

julia> g=Fun(CosSpace(),f.coefficients); # specify the coefficients directly

julia> f(cos(0.1))
2.70473560723178

julia> g(0.1)
2.7047356072317794
```

## Ultraspherical spaces

A key tool for solving differential equations are the ultraspherical spaces,
which can be defined by the span of derivatives of Chebyshev polynomials.  

Note that `Ultraspherical(1)` corresponds to the Chebyshev basis of the
second kind: $U_k(x) = {\sin (k+1) {\rm acos} x \over \sin {\rm acos} x}$.  
The relationship with Chebyshev polynomials follows from trigonemetric
identities: $T_k'(x) = k U_{k-1}(x)$.  

Converting between ultraspherical polynomials (with integer orders) is
extremely efficient: it requires $O(n)$ operations, where $n$ is the number of coefficients.

## Fourier and Laurent spaces

There are several different spaces to represent functions on periodic domains,
which are typically a `PeriodicInterval`, `Circle` or `PeriodicLine`.  

`CosSpace` represents expansion in cosine series:

$$f(\theta) = \sum_{k=0}^\infty f_k \cos k \theta$$

`SinSpace` represents expansion in sine series:

$$f(\theta) = \sum_{k=0}^\infty f_k \sin (k+1) \theta$$

`Taylor` represents expansion with only non-negative complex exponential terms:

$$f(\theta) = \sum_{k=0}^\infty f_k {\rm e}^{{\rm i} k \theta}$$

`Hardy{false}` represents expansion with only negative complex exponential terms:

$$f(\theta) = \sum_{k=0}^\infty f_k {\rm e}^{-{\rm i} (k+1) \theta}$$

`Fourier` represents functions that are sums of sines and cosines.  Note that
if a function has the form

$$f(\theta) = f_0 + \sum_{k=1}^\infty f_k^{\rm c} \cos k \theta + f_k^{\rm s} \sin k\theta$$

then the coefficients of the resulting `Fun` are order as $[f_0,f_1^{\rm s},f_1^{\rm c},…]$.
For example:

```jldoctest
julia> f = Fun(Fourier(),[1,2,3,4]);

julia> f(0.1)
4.979356652307978

julia> 1 + 2sin(0.1) + 3cos(0.1) + 4sin(2*0.1)
4.979356652307979
```

`Laurent` represents functions that are sums of complex exponentials.  Note that
if a function has the form

$$f(\theta) = \sum_{k=-\infty}^\infty f_k {\rm e}^{{\rm i} k \theta}$$

then the coefficients of the resulting `Fun` are order as $[f_0,f_{-1},f_1,…]$.
For example:

```jldoctest
julia> f = Fun(Laurent(),[1,2,3,4]);

julia> f(0.1)
9.895287137755096 - 0.694843906533417im

julia> 1 + 2exp(-im*0.1) + 3exp(im*0.1) + 4exp(-2im*0.1)
9.895287137755094 - 0.6948439065334167im
```



## Modifier spaces

Some spaces are built out of other spaces.  A simple example is `JacobiWeight(β,α,space)`
which weights `space`, which is typically `Chebyshev()` or `Jacobi(b,a)`,
 by a Jacobi weight `(1+x)^β*(1-x)^a`.


 ```@meta
 DocTestSetup = nothing
 ```


## Unset space

`UnsetSpace` is a special space that is used as a stand in when a
space has not yet been determined, particularly by operators.  
