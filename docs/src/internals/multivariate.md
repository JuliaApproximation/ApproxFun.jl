# Multivariate functions

Bivariate functions can also be represented in ApproxFun.  For example,

```julia
f=Fun((x,y)->exp(-x^2-y^2), ChebyshevInterval()^2)
```

constructs a function on the rectangle `[-1,1]^2`.  Just as in the 1D case, `f.coefficients` is a `Vector` of coefficients, and `f.space` is a `Space` that tells how to interpret the coefficients.

## TensorSpace

In the example above, the resulting space is a [`TensorSpace`](@ref).

## AbstractProductSpace

A `TensorSpace` is a subtype of `AbstractProductSpace`.  The purpose of `AbstractProductSpace`s beyond `TensorSpace` is that it is often convenient to use different bases in one-dimension depending on the order.  Thus we want to be able to represent functions in the basis `φ_k^j(x)ζ_j(y)`: for example, we could have the basis
```math
(1-x^2)^{m/2} \mathrm{P}_k^{(m,m)}(x) e^{imy}
```
(which is related to spherical harmonics), where ``\mathrm{P}_k^{(m,m)}(x)`` are Jacobi polynomials.

To handle this more general setting, an `AbstractProductSpace` implements `ApproxFunBase.columnspace`:  In the `(1-x^2)^(m/2) P_k^{(m,m)}(x) e^{imy}` example, we would have

```julia
columnspace(myproductspace,1)  ==  JacobiWeight(0,0,Legendre()),
columnspace(myproductspace,2)  ==  JacobiWeight(0.5,0.5,Jacobi(1,1)),
columnspace(myproductspace,k)  ==  JacobiWeight((k-1)/2,(k-1)/2,Jacobi(k-1,k-1))
```

The basis is then interlaced as in `TensorSpace`.

## ProductFun

A `Fun` is only one possible way to represent a bivariate function.  This mirrors linear algebra where a matrix can be represented in different formats depending on the usage `Matrix` or `SparseMatrixCSC`, for example.

Another format for representing bivariate functions in ApproxFun is a [`ProductFun`](@ref):  for example, the previous function could also be represented by

```julia
f=ProductFun((x,y)->exp(-x^2-y^2),ChebyshevInterval()^2)
```

A `ProductFun` also has two fields: `f.coefficients` and `f.space`, where `f.space` must be an `AbstractProductSpace`.   Here, `f.coefficients` is a `Vector{Fun{S,T}}` that represents a list of functions in `x`, where

```julia
space(f.coefficients[k]) == columnspace(f.space, k)
```

## LowRankFun

[`LowRankFun`](@ref) represents a function as a sum of outer products of 1D `Fun`s. This form tries to minimize the number of functions necessary by retaining the highest singular values of the function.
