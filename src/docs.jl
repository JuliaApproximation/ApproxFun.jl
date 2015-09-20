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
