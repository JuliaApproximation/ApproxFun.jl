# Library


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
Disk
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
PeriodicInterval
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
evaluate(::Space,::AbstractVector,::)
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
Jacobi
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
ApproxFun.values
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
chop
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
getindex(::Operator,::,::)
```

```@docs
getindex(::Operator,::Fun)
```


```@docs
\(::Operator,::)
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

## Bivariate

```@docs
LowRankFun
```
