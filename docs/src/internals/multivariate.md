Bivariate functions can also be represented in ApproxFun.  For example,

```julia
f=Fun((x,y)->exp(-x^2-y^2),Interval()^2)
```

constructs a function on the rectangle `[-1,1]^2`.  Just as in the 1D case, `f.coefficients` is a `Vector` of coefficients, and `f.space` is a `Space` that tells how to interpret the coefficients.

## TensorSpace

In the example above, the resulting space is a


## AbstractProductSpace

A `TensorSpace` is a subtype of `AbstractProductSpace`.  The purpose of `AbstractProductSpace`s beyond `TensorSpace` is that it is often convenient to use different bases in one-dimension depending on the order.  Thus we want to be able to represent functions in the basis `φ_k^j(x)ζ_j(y)`: for example, we could have the basis `(1-x^2)^j P_k^{(j,j)}(x) e^{ijy}` which is related to spherical harmonics, where `P_k^{(a,b)}` are the Jacobi polynomials.

To handle this more general setting, an `AbstractProductSpace` implements `columnspace`:  In the `(1-x^2)^j P_k^{(j,j)}(x) e^{ijy}` example, we would have

```julia
columnspace(myproductspace,1)  ==  JacobiWeight(0,0,Legendre()),
columnspace(myproductspace,2)  ==  JacobiWeight(1,1,Jacobi(1,1)),
columnspace(myproductspace,k)  ==  JacobiWeight(k-1,k-1,Jacobi(k-1,k-1))
```

The basis is then interlaced as in `TensorSpace`.

## ProductFun

A `Fun` is only one possible way to represent a bivariate function.  This mirrors linear algebra where a matrix can be represented in different formats depending on the usage `Matrix` or `SparseMatrixCSC`, for example.

Another format for representing bivariate functions in ApproxFun is a `ProductFun`:  for example, the previous function could also be represented by

```julia
f=ProductFun((x,y)->exp(-x^2-y^2),Interval()^2)
```

A `ProductFun` also has two fields: `f.coefficients` and `f.space`, where `f.space` must be an `AbstractProductSpace`.   Here, `f.coefficients` is a `Vector{Fun{S,T}}` that represents a list of functions in `x`, where

```julia
space(f.coefficients[k]) == columnspace(f.space,k)
```





## LowRankFun

`LowRankFun`

represents as a sum of outer products of 1D `Fun`s
