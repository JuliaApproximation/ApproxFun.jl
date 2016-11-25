# Constructors

`Fun`s in ApproxFun are instances of Julia types with one field to store coefficients and another
to describe the function space. Similarly, each function space has one field describing
its domain, or another function space. Let's explore:

```julia
x = Fun(identity,-1..1)
f = exp(x)
g = f/sqrt(1-x^2)
space(f)   # Chebyshev(Interval(-1.0,1.0))
space(g)   # JacobiWeight(-0.5,-0.5,Interval(-1.0,1.0))
```

The absolute value is
another case where the space of the output is inferred from the operation:

```julia
f = Fun(x->cospi(5x),-1..1)
g = abs(f)
space(f)   # Chebyshev(Interval(-1.0,1.0))
space(g)   # PiecewiseSpace((Chebyshev(Interval(-1.0,-0.9)),...))
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
```julia
f = Fun(Taylor(),[1,2,3])  # represents 1 + 2z + 3z^2
```



## Using ApproxFun for “manual” interpolation

The ApproxFun package for Julia implements all of the necessary operations for Chebyshev interpolation and operations (like differentiation or integration) on Chebyshev interpolants.

Normally, you give it a function f and a domain d, and construct the Chebyshev interpolant by `fc = Fun(f, d)`. The ApproxFun package figures out the necessary number of Chebyshev points (i.e., the polynomial order) required to interpolate f to nearly machine precision, so that subsequent operations on fc can be viewed as "exact".

However, in cases where the function to be interpolated is extremely expensive, and possibly even is evaluated by an external program, it is convenient to be able to decide on the desired Chebyshev order in advance, evaluate the function at those points "manually", and then construct the Chebyshev interpolant. However, this procedure isn't documented in the ApproxFun manual. In this notebook, we show how to do that for the example `f(x) = exp(2x)` on the domain `0..1`.

### Manually constructing the Chebyshev interpolant

```julia
using ApproxFun
f(x) = exp(2x) # the function to be interpolated
d = Space(0..1) # the domain to interpolate in
```

Let's interpolate using only `N = 10` points. We'll start by using the ApproxFun.points function to get the x coordinates to interpolate at:

```julia
N = 10
x = points(d, N)
10-element Array{Float64,1}:
 0.00615583
 0.0544967
 0.146447
 0.273005
 0.421783
 0.578217
 0.726995
 0.853553
 0.945503
 0.993844
```

Note that these are Chebyshev points, clustered near the endpoints of the domain.

Next, we'll evaluate our function f at these points:
```julia
y = f(x)
10-element Array{Float64,1}:
 1.01239
 1.11516
 1.3403
 1.72635
 2.32464
 3.17858
 4.28016
 5.51299
 6.62603
 7.29864
```

Finally, we'll construct the Chebyshev interpolant fc by using the lowest-level Fun constructor, where we pass in an array of Chebyshev coefficients along with the domain. To get the Chebyshev coefficients, we'll use the ApproxFun.transform method:

```julia
fc = Fun(ApproxFun.transform(d, y), d)
Fun([3.441523869125335,3.0725234451419356,0.7380008479667985,0.12052005327473987,0.01488052831835938,0.0014758267278678237,0.00012226103967574454,8.694251606738979e-6,5.415128415009462e-7,2.993315462163082e-8],Chebyshev(【0.0,1.0】))
```

### Using the Chebyshev interpolant

We can evaluate the interpolation at a point x ∈ d just by calling fc(x). Since we only used 10 points, it is not accurate to machine precision, but it is still pretty accurate:

```julia
fc(0.2) - f(0.2)
1.3975267609822595e-9
```

We can compute the derivative at an arbitrary point simply by fc'(x), which will give us a good approximation for the exact derivative 2*exp(2x):

```julia
fc'(0.2) - 2*exp(2*0.2)
5.5879429972094385e-9
```
