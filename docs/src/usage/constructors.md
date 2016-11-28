# Constructors

`Fun`s in ApproxFun are instances of Julia types with one field to store coefficients and another
to describe the function space. Similarly, each function space has one field describing
its domain, or another function space. Let's explore:


```@meta
DocTestSetup = quote
    using ApproxFun
end
```


```jldoctest
julia> x = Fun(identity,-1..1);

julia> f = exp(x);

julia> g = f/sqrt(1-x^2);

julia> space(f)   # Output is pretty version of Chebyshev(Interval(-1.0,1.0))
Chebyshev(【-1.0,1.0】)

julia> space(g)   # Output is pretty version of  JacobiWeight(-0.5,-0.5,Interval(-1.0,1.0))
(1-x^2)^-0.5[Chebyshev(【-1.0,1.0】)]
```

The absolute value is
another case where the space of the output is inferred from the operation:

```jldoctest
julia> f = Fun(x->cospi(5x),-1..1);

julia> g = abs(f);

julia> space(f)   
Chebyshev(【-1.0,1.0】)

julia> space(g)
Chebyshev(【-1.0,-0.9000000000000002】)⨄Chebyshev(【-0.9000000000000002,-0.6999999999999996】)⨄Chebyshev(【-0.6999999999999996,-0.5000000000000001】)⨄Chebyshev(【-0.5000000000000001,-0.30000000000000043】)⨄Chebyshev(【-0.30000000000000043,-0.09999999999999962】)⨄Chebyshev(【-0.09999999999999962,0.10000000000000053】)⨄Chebyshev(【0.10000000000000053,0.29999999999999966】)⨄Chebyshev(【0.29999999999999966,0.500000000000001】)⨄Chebyshev(【0.500000000000001,0.6999999999999998】)⨄Chebyshev(【0.6999999999999998,0.9000000000000006】)⨄Chebyshev(【0.9000000000000006,1.0】)
```

## Convenience constructors

The default space is `Chebyshev`, which can represent non-periodic functions on intervals.  Each `Space` type has a default domain: for `Chebyshev` this is `-1..1`, for Fourier and Laurent this is `-π..π`.  Thus the following
are synonyms:
```julia
Fun(exp,Chebyshev(Interval(-1,1)))
Fun(exp,Chebyshev(Interval()))
Fun(exp,Chebyshev(-1..1))
Fun(exp,Chebyshev())
Fun(exp,-1..1)
Fun(exp,Interval())
Fun(exp,Interval(-1,1))
Fun(exp)
```
If a function is not specified, then it is taken to be `identity`.  Thus we have the
following synonyms:
```julia
x = Fun(identity,-1..1)
x = Fun(-1..1)
x = Fun(identity)
x = Fun()
```


## Specifying coefficients explicitly

It is sometimes necessary to specify coefficients explicitly.  This is possible
via specifying the space followed by a vector of coefficients:
```jldoctest
julia> f = Fun(Taylor(),[1,2,3]);  # represents 1 + 2z + 3z^2

julia> f(0.1)
1.23

julia> 1+2*0.1+3*0.1^2
1.23
```



## Using ApproxFun for “manual” interpolation

The ApproxFun package for Julia implements all of the necessary operations for Chebyshev interpolation and operations (like differentiation or integration) on Chebyshev interpolants.

Normally, you give it a function f and a domain d, and construct the Chebyshev interpolant by `fc = Fun(f, d)`. The ApproxFun package figures out the necessary number of Chebyshev points (i.e., the polynomial order) required to interpolate f to nearly machine precision, so that subsequent operations on fc can be viewed as "exact".

However, in cases where the function to be interpolated is extremely expensive, and possibly even is evaluated by an external program, it is convenient to be able to decide on the desired Chebyshev order in advance, evaluate the function at those points "manually", and then construct the Chebyshev interpolant. However, this procedure isn't documented in the ApproxFun manual. In this notebook, we show how to do that for the example `f(x) = exp(2x)` on the domain `0..1`.
