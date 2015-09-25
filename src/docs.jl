



## Fun.jl docs

doc"""
    domain(::Fun)

returns the domain that a `Fun` is defined on
"""
domain(fun::Fun)


doc"""
    setdomain(fun,domain)

returns `fun` projected onto `domain`
"""
domain(fun::Fun,domain::Domain)


doc"""
    space(fun)

returns the space of `fun`
"""
space(fun::Fun)



doc"""
    values(fun)

returns `fun` evaluated at `points(fun)`
"""
Base.values(fun::Fun)




doc"""
    points(fun)

returns a grid of points that the fun can be transformed into values
and back
"""
points(fun::Fun)


doc"""
    length(fun) -> Integer

returns the number of coefficients of a fun
"""
Base.length(fun::Fun)

doc"""
    stride(fun)

returns the stride of the coefficients, checked
numerically
"""
Base.stride(fun::Fun)



doc"""
    reverseorientation(fun)

return `fun` on a reversed orientated contour
"""
reverseorientation(f::Fun)



## Operator.jl docs

doc"""
    `Operator{T}` represents a general infinite operator
"""
Operator

doc"""
    `Functional{T}` represents a row operator
"""
Functional

doc"""
    `InfiniteOperator{T}` represents an operator with an infinite number of rows
"""
InfiniteOperator

doc"""
    `BandedBelowOperator{T}` represents an operator banded banded below.
    The bandwidth can be found with bandinds(op,1).
"""
BandedBelowOperator

doc"""
    `AlmostBandedOperator{T}` represents an operator that is banded apart from
    a finite number of dense rows. The bandwidth can be found with bandinds(op).
"""
AlmostBandedOperator

doc"""
    `BandedOperator{T}` represents a banded operator. The bandwidth can be found with bandinds(op).
"""
BandedOperator
