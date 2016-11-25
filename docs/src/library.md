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
Interval
```

```@docs
Line
```

```@docs
PeriodicInterval
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
evaluate(::Space,::Vector,::)
```

```@docs
ApproxFun.dimension(::Space)
```


## Included spaces

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
Taylor
```

```@docs
Hardy{false}
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
Jacobi
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
domainspace
```

```@docs
rangespace
```

```@docs
Base.getindex(::Operator,::,::)
```

```@docs
linsolve
```

```@docs
Conversion
```

```@docs
qrfact(::Operator)
```

```@docs
cache(::Operator)
```
