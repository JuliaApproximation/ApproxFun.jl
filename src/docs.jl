



## Fun.jl docs

# Constructors

doc"""
    Fun(s::Space,coefficients::AbstractVector)

returns a `Fun` with the specified `coefficients` in the space `s`
"""
Fun(::Space,::AbstractVector)

doc"""
    Fun(f,s::Space)

return a `Fun` representing the function, number, or vector `f` in the
space `s`.  If `f` is vector-valued, it returns a vector-valued analogue
of `s`.
"""
Fun(_,::Space)

doc"""
    Fun(f,d::Domain)

returns `Fun(f,Space(d))`, that is, it uses the default space for the specified
domain.
"""
Fun(_,::Domain)


doc"""
    Fun(s::Space)

returns `Fun(identity,s)`
"""
Fun(::Space)

doc"""
    Fun(f)

returns `Fun(f,Chebyshev())`
"""
Fun(f)

doc"""
    Fun()

returns `Fun(identity,Chebyshev())`.
"""
Fun()

doc"""
    ones(d::Space)

Return the `Fun` that represents the function one on the specified space.
"""
Base.ones(::Space)

doc"""
    zeros(d::Space)

Return the `Fun` that represents the function one on the specified space.
"""
Base.zeros(::Space)

# accessors

doc"""
    domain(f::Fun)

returns the domain that `f` is defined on
"""
domain(fun::Fun)


doc"""
    setdomain(f::Fun,d::Domain)

returns `f` projected onto `domain`
"""
setdomain(::Fun,::Domain)


doc"""
    space(f::Fun)

returns the space of `f`
"""
space(f::Fun)



doc"""
    values(f::Fun)

returns `f` evaluated at `points(f)`
"""
Base.values(::Fun)




doc"""
    points(f::Fun)

returns a grid of points that `f` can be transformed into values
and back
"""
points(::Fun)

doc"""
    points(s::Space,n::Integer)

returns a grid of approximately `n` points, for which a transform exists
from values at the grid to coefficients in the space `s`.
"""
points(::Space,::Integer)

doc"""
    extrapolate(f::Fun,x)

returns an extrapolation of `f` from its domain to `x`.
"""
extrapolate(::Fun,x)


doc"""
    coefficients(f::Fun) -> Vector

returns the coefficients of `f`, corresponding to the space `space(f)`.
"""
coefficients(::Fun)


"""
    coefficients(f::Fun,s::Space) -> Vector

returns the coefficients of `f` in the space `s`, which
may not be the same as `space(f)`.
"""
coefficients(::Fun,::Space)

"""
    coefficients(cfs::AbstractVector,fromspace::Space,tospace::Space) -> Vector

converts coefficients in `fromspace` to coefficients in `tospace`
"""
coefficients(::AbstractVector,::Space,::Space)


doc"""
    ncoefficients(f::Fun) -> Integer

returns the number of coefficients of a fun
"""
ncoefficients(::Fun)

doc"""
    stride(f::Fun)

returns the stride of the coefficients, checked
numerically
"""
Base.stride(::Fun)



## Modifiers

doc"""
   chop(f::Fun,tol) -> Fun

reduces the number of coefficients by dropping the tail that is below the specified tolerance.
"""
chop(::Fun,_)

doc"""
    reverseorientation(f::Fun)

return `f` on a reversed orientated contour.
"""
reverseorientation(::Fun)


## Spaces

doc"""
    canonicalspace(s::Space)

returns a space that is used as a default to implement missing functionality,
e.g., evaluation.  Implement a `Conversion` operator or override `coefficients` to support this.
"""
ApproxFun.canonicalspace(::Space)

doc"""
    transform(s::Space,vals::Vector)

Transform values on the grid specified by `points(s,length(vals))` to coefficients in the space `s`.
Defaults to `coefficients(transform(canonicalspace(space),values),canonicalspace(space),space)`
"""
transform(::Space,::Vector)

doc"""
    itransform(s::Space,coefficients::AbstractVector)

Transform coefficients back to values.  Defaults to using `canonicalspace` as in `transform`.
"""
itransform(::Space, ::AbstractVector)


doc"""
    evaluate(coefficients::AbstractVector, sp::Space, x)

Evaluates the expansion at a point `x`.
If `x` is in the domain, then this should return zero.
"""
evaluate(::AbstractVector, ::Space, _)



doc"""
    spacescompatible

Specifies equality of spaces while also supporting `AnyDomain`.
"""
spacescompatible(::Space,::Space)

doc"""
    conversion_type(a::Space,b::Space)

returns a `Space` that has a banded conversion operator to both `a` and `b`.
Override `ApproxFun.conversion_rule` when adding new `Conversion` operators.
"""
conversion_type(::Space,::Space)

doc"""
    dimension(s::Space)

returns the dimension of `s`, which is the maximum number of coefficients.
"""
dimension(::Space)


## Operator.jl docs

doc"""
    Operator{T}

is an abstract type to represent linear operators between spaces.
"""
Operator

doc"""
    domainspace(op::Operator)

gives the domain space of `op`.  That is, `op*f` will first convert `f` to
a `Fun` in the space `domainspace(op)` before applying the operator.
"""
domainspace(::Operator)

doc"""
    rangespace(op::Operator)

gives the range space of `op`.  That is, `op*f` will return a `Fun` in the
space `rangespace(op)`, provided `f` can be converted to a `Fun` in
`domainspace(op)`.
"""
rangespace(::Operator)


doc"""
    bandwidths(op::Operator)

returns the bandwidth of `op` in the form `(l,u)`, where `l ≥ 0` represents
the number of subdiagonals and `u ≥ 0` represents the number of superdiagonals.
"""
bandwidths(::Operator)

doc"""
    promotedomainspace(S::Operator,sp::Space)

returns the operator `S` but acting on the space `sp`.
"""
promotedomainspace(::Operator,::Space)

doc"""
    promoterangespace(S::Operator,sp::Space)

returns the operator `S` acting on the same space, but now return
functions in the specified range space `sp`
"""
promoterangespace(::Operator,::Space)

doc"""
    choosedomainspace(S::Operator,rangespace::Space)

returns a space `ret` so that `promotedomainspace(S,ret)` has the
specified range space.
"""
choosedomainspace(::Operator,::Space)


doc"""
    op[k,j]

returns the `k`th coefficient of `op*Fun([zeros(j-1);1],domainspace(op))`.
"""
Base.getindex(::Operator,k,j)


doc"""
    op[f::Fun]

constructs the operator `op*Multiplication(f)`, that is, it multiplies on the right
by `f` first.  Note that `op*f` is different: it applies `op` to `f`.
"""
Base.getindex(::Operator,::Fun)


doc"""
    Conversion(fromspace::Space,tospace::Space)

represents a conversion operator between `fromspace` and `tospace`, when available.
"""
Conversion(::Space,::Space)
