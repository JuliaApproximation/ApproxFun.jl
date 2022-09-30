```@meta
DocTestSetup  = quote
    using ApproxFun, LinearAlgebra
end
```

# Constructors

`Fun`s in ApproxFun are instances of Julia types with one field to store coefficients and another to describe the function space. Similarly, each function space has one field describing its domain, or another function space. Let's explore:

```@repl
x = Fun(identity,-1..1);
f = exp(x);
g = f/sqrt(1-x^2);
space(f)  # Output is pretty version of Chebyshev(Interval(-1.0,1.0))
space(g)  # Output is pretty version of JacobiWeight(-0.5,-0.5,Interval(-1.0,1.0))
```

The absolute value is another case where the space of the output is inferred from the operation:

```@repl
f = Fun(x->cospi(5x),-1..1);
g = abs(f);
space(f)
space(g)
```

## Convenience constructors

The default space is `Chebyshev`, which can represent non-periodic functions on intervals.  Each `Space` type has a default domain: for `Chebyshev` this is `-1..1`, for Fourier and Laurent this is `-π..π`.  Thus the following are equivalent:

```julia
Fun(exp,Chebyshev(Interval(-1,1)))
Fun(exp,Chebyshev(ChebyshevInterval()))
Fun(exp,Chebyshev(-1..1))
Fun(exp,Chebyshev())
Fun(exp,-1..1)
Fun(exp,ChebyshevInterval())
Fun(exp,Interval(-1,1))
Fun(exp)
```

If a function is not specified, then it is taken to be `identity`.  Thus we have the following equivalent constructions:

```julia
x = Fun(identity, -1..1)
x = Fun(-1..1)
x = Fun(identity)
x = Fun()
```

## Specifying coefficients explicitly

It is sometimes necessary to specify coefficients explicitly.  This is possible via specifying the space followed by a vector of coefficients:

```@repl
f = Fun(Taylor(), [1,2,3]);  # represents 1 + 2z + 3z^2
f(0.1)
1 + 2*0.1 + 3*0.1^2
```

In higher dimensions, ApproxFun will sum products of the 1D basis functions. So if ``\mathop{T}_i(x)`` is the ``i``th basis function, then a 2D function can be approximated as the following:

```math
\mathop{f}(x,y) = \sum_{i,j} c_{ij} \mathop{T}_i(x) \mathop{T}_j(y).
```

The products will be ordered lexicographically by the degree of the polynomial, i.e., in the order

```math
\mathop{T}_0(x) \mathop{T}_0(y),\ \mathop{T}_0(x) \mathop{T}_1(y),\ \mathop{T}_1(x) \mathop{T}_0(y),\ \mathop{T}_0(x) \mathop{T}_2(y),\ \mathop{T}_1(x) \mathop{T}_1(y),\ \mathop{T}_2(x) \mathop{T}_0(y),\ ….
```

For example, if we are in the two dimensional CosSpace space and we have coefficients ``\{c_1, c_2, c_3\}``, then

```math
\mathop{f}(x, y) = c_1 \cos(0 x) \cos(0 y) + c_2 \cos(0 x) \cos(1 y) + c_3 \cos(1 x) \cos(0 y).
```

This is illustrated in the following code:

```@repl
f = Fun(CosSpace()^2, [1,2,3])
f(1,2)
1cos(0*1)*cos(0*2) + 2cos(0*1)*cos(1*2) + 3cos(1*1)*cos(0*2)
```

## Using ApproxFun for “manual” interpolation

The ApproxFun package for Julia implements all of the necessary operations for Chebyshev interpolation and operations (like differentiation or integration) on Chebyshev interpolants.

Normally, you give it a function f and a domain d, and construct the Chebyshev interpolant by `fc = Fun(f, d)`. The ApproxFun package figures out the necessary number of Chebyshev points (i.e., the polynomial order) required to interpolate f to nearly machine precision, so that subsequent operations on fc can be viewed as "exact".

However, in cases where the function to be interpolated is extremely expensive, and possibly even is evaluated by an external program, it is convenient to be able to decide on the desired Chebyshev order in advance, evaluate the function at those points "manually", and then construct the Chebyshev interpolant. An example showing how to do this is given in the [ApproxFun FAQ](../faq.md).
