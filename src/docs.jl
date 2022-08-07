



## Fun.jl docs

# Constructors

"""
    Fun(s::Space, coefficients::AbstractVector)

Return a `Fun` with the specified `coefficients` in the space `s`

# Examples
```jldoctest
julia> f = Fun(Fourier(), [1,1]);

julia> f(0.1) == 1 + sin(0.1)
true

julia> f = Fun(Chebyshev(), [1,1]);

julia> f(0.1) == 1 + 0.1
true
```
"""
Fun(::Space,::AbstractVector)

"""
    Fun(f, s::Space)

Return a `Fun` representing the function, number, or vector `f` in the
space `s`.  If `f` is vector-valued, it Return a vector-valued analogue
of `s`.

# Examples
```jldoctest
julia> f = Fun(x->x^2, Chebyshev())
Fun(Chebyshev(), [0.5, 0.0, 0.5])

julia> f(0.1) == (0.1)^2
true
```
"""
Fun(_, ::Space)

"""
    Fun(f, d::Domain)

Return `Fun(f, Space(d))`, that is, it uses the default space for the specified
domain.

# Examples
```jldoctest
julia> f = Fun(x->x^2, 0..1)
Fun(Chebyshev(0..1), [0.375, 0.5, 0.12499999999999997])

julia> f(0.1) ≈ (0.1)^2
true
```
"""
Fun(_, ::Domain)


"""
    Fun(s::Space)

Return `Fun(identity,s)`

# Examples
```jldoctest
julia> x = Fun(Chebyshev())
Fun(Chebyshev(), [0.0, 1.0])

julia> x(0.1)
0.1
```
"""
Fun(::Space)

"""
    Fun(f)

Return `Fun(f, Chebyshev())`

# Examples
```jldoctest
julia> f = Fun(x->x^2)
Fun(Chebyshev(), [0.5, 0.0, 0.5])

julia> f(0.1) == (0.1)^2
true
```
"""
Fun(f)

"""
    Fun()

Return `Fun(identity, Chebyshev())`, which represents the identity function in `-1..1`.

# Examples
```jldoctest
julia> f = Fun(Chebyshev())
Fun(Chebyshev(), [0.0, 1.0])

julia> f(0.1)
0.1
```
"""
Fun()

"""
    ones(d::Space)

Return the `Fun` that represents the function one on the specified space.

# Examples
```jldoctest
julia> ones(Chebyshev())
Fun(Chebyshev(), [1.0])
```
"""
ones(::Space)

"""
    zeros(d::Space)

Return the `Fun` that represents the function one on the specified space.

# Examples
```jldoctest
julia> zeros(Chebyshev())
Fun(Chebyshev(), [0.0])
```
"""
zeros(::Space)

# accessors

"""
    domain(f::Fun)

Return the domain that `f` is defined on.

# Examples
```jldoctest
julia> f = Fun(x->x^2)
Fun(Chebyshev(), [0.5, 0.0, 0.5])

julia> domain(f)
-1.0..1.0 (Chebyshev)
```
"""
domain(fun::Fun)


"""
    setdomain(f::Fun, d::Domain)

Return `f` projected onto `domain`.

!!! note
    The new function may differ from the original one, as the coefficients are left unchanged.

# Examples
```jldoctest
julia> f = Fun(x->x^2)
Fun(Chebyshev(), [0.5, 0.0, 0.5])

julia> setdomain(f, 0..1)
Fun(Chebyshev(0..1), [0.5, 0.0, 0.5])
```
"""
setdomain(::Fun,::Domain)


"""
    space(f::Fun)

Return the space of `f`.

# Examples
```jldoctest
julia> f = Fun(x->x^2)
Fun(Chebyshev(), [0.5, 0.0, 0.5])

julia> space(f)
Chebyshev()
```
"""
space(f::Fun)



"""
    values(f::Fun)

Return `f` evaluated at `points(f)`.

# Examples
```jldoctest
julia> f = Fun(x->x^2)
Fun(Chebyshev(), [0.5, 0.0, 0.5])

julia> values(f)
3-element Vector{Float64}:
 0.75
 0.0
 0.75

julia> map(x->x^2, points(f)) ≈ values(f)
true
```
"""
values(::Fun)




"""
    points(f::Fun)

Return a grid of points that `f` can be transformed into values
and back.

# Examples
```jldoctest
julia> f = Fun(x->x^2)
Fun(Chebyshev(), [0.5, 0.0, 0.5])

julia> points(f)
3-element Vector{Float64}:
  0.8660254037844386
  0.0
 -0.8660254037844386
```
"""
points(::Fun)

"""
    points(s::Space,n::Integer)

Return a grid of approximately `n` points, for which a transform exists
from values at the grid to coefficients in the space `s`.

# Examples
```jldoctest
julia> points(Chebyshev(), 4)
4-element Vector{Float64}:
  0.9238795325112867
  0.3826834323650898
 -0.3826834323650898
 -0.9238795325112867
```
"""
points(::Space,::Integer)

"""
    extrapolate(f::Fun,x)

Return an extrapolation of `f` from its domain to `x`.

# Examples
```jldoctest
julia> f = Fun(x->x^2)
Fun(Chebyshev(), [0.5, 0.0, 0.5])

julia> domain(f)
-1.0..1.0 (Chebyshev)

julia> extrapolate(f, 2)
4.0
```
"""
extrapolate(::Fun,x)


"""
    coefficients(f::Fun) -> Vector

Return the coefficients of `f`, corresponding to the space `space(f)`.

# Examples
```jldoctest
julia> f = Fun(x->x^2)
Fun(Chebyshev(), [0.5, 0.0, 0.5])

julia> coefficients(f)
3-element Vector{Float64}:
 0.5
 0.0
 0.5
```
"""
coefficients(::Fun)


"""
    coefficients(f::Fun, s::Space) -> Vector

Return the coefficients of `f` in the space `s`, which
may not be the same as `space(f)`.

# Examples
```jldoctest
julia> f = Fun(x->x^2)
Fun(Chebyshev(), [0.5, 0.0, 0.5])

julia> coefficients(f, NormalizedChebyshev())
3-element Vector{Float64}:
 0.8862269254527579
 0.0
 0.6266570686577501

julia> coefficients(f, Legendre())
3-element Vector{Float64}:
 0.33333333333333337
 0.0
 0.6666666666666666
```
"""
coefficients(::Fun,::Space)

"""
    coefficients(cfs::AbstractVector, fromspace::Space, tospace::Space) -> Vector

Convert coefficients in `fromspace` to coefficients in `tospace`

# Examples
```jldoctest
julia> f = Fun(x->x^2)
Fun(Chebyshev(), [0.5, 0.0, 0.5])

julia> coefficients(f, Chebyshev(), Legendre())
3-element Vector{Float64}:
 0.33333333333333337
 0.0
 0.6666666666666666

julia> Fun(x->x^2, Legendre())
Fun(Legendre(), [0.33333333333333337, 0.0, 0.6666666666666666])
```
"""
coefficients(::AbstractVector,::Space,::Space)


"""
    ncoefficients(f::Fun) -> Integer

Return the number of coefficients of a fun

# Examples
```jldoctest
julia> f = Fun(x->x^2)
Fun(Chebyshev(), [0.5, 0.0, 0.5])

julia> ncoefficients(f)
3
```
"""
ncoefficients(::Fun)

"""
    stride(f::Fun)

Return the stride of the coefficients, checked numerically
"""
stride(::Fun)



## Modifiers

"""
    chop(f::Fun, tol) -> Fun

Reduce the number of coefficients by dropping the tail that is below the specified tolerance.

# Examples
```jldoctest
julia> f = Fun(Chebyshev(), [1,2,3,0,0,0])
Fun(Chebyshev(), [1, 2, 3, 0, 0, 0])

julia> chop(f)
Fun(Chebyshev(), [1, 2, 3])
```
"""
chop(::Fun,_)

"""
    reverseorientation(f::Fun)

Return `f` on a reversed orientated contour.
"""
reverseorientation(::Fun)


## Spaces

"""
    canonicalspace(s::Space)

Return a space that is used as a default to implement missing functionality,
e.g., evaluation.  Implement a `Conversion` operator or override `coefficients` to support this.

# Examples
```jldoctest
julia> f = Fun(x->x^2, NormalizedLegendre())
Fun(NormalizedLegendre(), [0.4714045207910318, 0.0, 0.4216370213557839])

julia> ApproxFunBase.canonicalspace(f)
Legendre()
```
"""
ApproxFun.canonicalspace(::Space)

"""
    transform(s::Space, vals::Vector)

Transform values on the grid specified by `points(s,length(vals))` to coefficients in the space `s`.
Defaults to `coefficients(transform(canonicalspace(space),values),canonicalspace(space),space)`

# Examples
```jldoctest
julia> v = map(x -> x^2, points(Chebyshev(), 4));

julia> transform(Chebyshev(), v)
4-element Vector{Float64}:
 0.5
 0.0
 0.5
 0.0
```
"""
transform(::Space,::Vector)

"""
    itransform(s::Space,coefficients::AbstractVector)

Transform coefficients back to values.  Defaults to using `canonicalspace` as in `transform`.

# Examples
```jldoctest
julia> v = itransform(Chebyshev(), [0.5, 0, 0.5])
3-element Vector{Float64}:
 0.75
 0.0
 0.75

julia> f = Fun(x->x^2, Chebyshev())
Fun(Chebyshev(), [0.5, 0.0, 0.5])

julia> values(f)
3-element Vector{Float64}:
 0.75
 0.0
 0.75
```
"""
itransform(::Space, ::AbstractVector)


"""
    evaluate(coefficients::AbstractVector, sp::Space, x)

Evaluates the expansion at a point `x`.
If `x` is in the domain, then this should return zero.
"""
evaluate(::AbstractVector, ::Space, _)



"""
    spacescompatible

Specifies equality of spaces while also supporting `AnyDomain`.
"""
spacescompatible(::Space,::Space)

"""
    conversion_type(a::Space,b::Space)

Return a `Space` that has a banded conversion operator to both `a` and `b`.
Override `ApproxFun.conversion_rule` when adding new `Conversion` operators.
"""
conversion_type(::Space,::Space)

"""
    dimension(s::Space)

Return the dimension of `s`, which is the maximum number of coefficients.
"""
dimension(::Space)


## Operator.jl docs

"""
    Operator{T}

is an abstract type to represent linear operators between spaces.
"""
Operator

"""
    domainspace(op::Operator)

gives the domain space of `op`.  That is, `op*f` will first convert `f` to
a `Fun` in the space `domainspace(op)` before applying the operator.
"""
domainspace(::Operator)

"""
    rangespace(op::Operator)

gives the range space of `op`.  That is, `op*f` will return a `Fun` in the
space `rangespace(op)`, provided `f` can be converted to a `Fun` in
`domainspace(op)`.
"""
rangespace(::Operator)


"""
    bandwidths(op::Operator)

Return the bandwidth of `op` in the form `(l,u)`, where `l ≥ 0` represents
the number of subdiagonals and `u ≥ 0` represents the number of superdiagonals.
"""
bandwidths(::Operator)

"""
    promotedomainspace(S::Operator,sp::Space)

Return the operator `S` but acting on the space `sp`.
"""
promotedomainspace(::Operator,::Space)

"""
    promoterangespace(S::Operator,sp::Space)

Return the operator `S` acting on the same space, but now return
functions in the specified range space `sp`
"""
promoterangespace(::Operator,::Space)

"""
    choosedomainspace(S::Operator,rangespace::Space)

Return a space `ret` so that `promotedomainspace(S,ret)` has the
specified range space.
"""
choosedomainspace(::Operator,::Space)


"""
    op[k,j]

Return the `k`th coefficient of `op*Fun([zeros(j-1);1],domainspace(op))`.
"""
getindex(::Operator,k,j)


"""
    op[f::Fun]

constructs the operator `op * Multiplication(f)`, that is, it multiplies on the right
by `f` first.  Note that `op * f` is different: it applies `op` to `f`.

# Examples
```jldoctest
julia> x = Fun()
Fun(Chebyshev(), [0.0, 1.0])

julia> D = Derivative()
ConcreteDerivative : ApproxFunBase.UnsetSpace() → ApproxFunBase.UnsetSpace()

julia> D2 = D[x]
TimesOperator : ApproxFunBase.UnsetSpace() → ApproxFunBase.UnsetSpace()

julia> twox = D2 * x
Fun(Ultraspherical(1), [0.0, 1.0])

julia> twox(0.1) ≈ 2 * 0.1
true
```
"""
getindex(::Operator,::Fun)


"""
    Conversion(fromspace::Space,tospace::Space)

Represent a conversion operator between `fromspace` and `tospace`, when available.
"""
Conversion(::Space,::Space)
