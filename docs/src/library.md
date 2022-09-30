# Library

## Domains
```@docs
Arc
```

```@docs
Circle
```

```@docs
Curve
```

```@docs
ApproxFunFourier.Disk
```

```@docs
Segment
```

```@docs
Interval
```

```@docs
Line
```

```@docs
PeriodicSegment
```

```@docs
ApproxFun.Point
```

```@docs
ProductDomain
```

```@docs
Ray
```

```@docs
UnionDomain
```

```@docs
∂
```

## Accessing information about a spaces

```@docs
ApproxFun.canonicalspace
```

```@docs
itransform
```

```@docs
transform
```

```@docs
evaluate
```

```@docs
ApproxFun.dimension(::Space)
```


## Inbuilt spaces

```@docs
SequenceSpace
```

```@docs
ConstantSpace
```

```@docs
Chebyshev
```

```@docs
Hermite
```

```@docs
Jacobi
```

```@docs
Laguerre
```

```@docs
Ultraspherical
```

```@docs
Taylor
```

```@docs
Hardy
```

```@docs
Fourier
```

```@docs
Laurent
```

```@docs
CosSpace
```

```@docs
SinSpace
```

```@docs
JacobiWeight
```

```@docs
ApproxFun.LogWeight
```

```@docs
ApproxFun.ArraySpace
```

```@docs
TensorSpace
```

## Constructing a Fun

```@docs
Fun
```

```@docs
ones(::Space)
```

```@docs
zeros(::Space)
```

## Accessing information about a Fun

```@docs
domain
```

```@docs
coefficients
```

```@docs
extrapolate
```

```@docs
ncoefficients
```

```@docs
points
```

```@docs
space
```

```@docs
values(::Fun)
```

```@docs
stride(::Fun)
```

## Modify a Fun


```@docs
reverseorientation
```

```@docs
ApproxFun.setdomain
```

```@docs
chop(::Fun, ::Any)
```

## Bivariate Fun

```@docs
LowRankFun
```


## Operators

```@docs
Operator
```

```@docs
BandedMatrices.bandwidths(::Operator)
```

```@docs
domainspace
```

```@docs
rangespace
```

```@docs
getindex(::Operator,::Any,::Any)
```

```@docs
getindex(::Operator,::Fun)
```


```@docs
\(::Operator,::Any)
```

```@docs
qr(::Operator)
```

```@docs
cache(::Operator)
```



## Inbuilt operators

```@docs
Conversion
```

```@docs
Derivative
```

```@docs
Dirichlet
```

```@docs
Evaluation
```

```@docs
Integral
```

```@docs
Laplacian
```

```@docs
Multiplication
```

```@docs
Neumann
```

```@docs
PartialInverseOperator
```
