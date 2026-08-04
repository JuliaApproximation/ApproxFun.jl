module ApproxFunVectorInterfaceExt

using ApproxFun: Fun, space, coefficients, innerproduct
using ApproxFun.ApproxFunBase: spacescompatible
using VectorInterface: VectorInterface, One, scalartype, zerovector, zerovector!,
    scale, scale!, scale!!, add, add!

# Funs in compatible spaces share a coefficient basis, so the vector operations
# can act directly on the (possibly differently truncated) coefficient vectors.
# The in-place (`!`) variants require compatible spaces and mutable coefficient
# vectors; the `!!` variants fall back to the out-of-place versions otherwise.

VectorInterface.scalartype(::Type{<:Fun{<:Any,T}}) where {T} = scalartype(T)

function _check_compatible(y::Fun, x::Fun)
    spacescompatible(y, x) || throw(ArgumentError(
        "in-place operations require Funs in compatible spaces, got $(space(y)) and $(space(x))"))
    return nothing
end

_resizable(c) = c isa Vector

function VectorInterface.zerovector(f::Fun, ::Type{S}) where {S<:Number}
    return Fun(space(f), zerovector(coefficients(f), S))
end
VectorInterface.zerovector!(f::Fun) = (zerovector!(coefficients(f)); f)
function VectorInterface.zerovector!!(f::Fun)
    c = coefficients(f)
    c′ = VectorInterface.zerovector!!(c)
    return c′ === c ? f : Fun(space(f), c′)
end

VectorInterface.scale(f::Fun, α::Number) = Fun(space(f), scale(coefficients(f), α))
VectorInterface.scale!(f::Fun, α::Number) = (scale!(coefficients(f), α); f)
function VectorInterface.scale!!(f::Fun, α::Number)
    c = coefficients(f)
    c′ = scale!!(c, α)
    return c′ === c ? f : Fun(space(f), c′)
end

function VectorInterface.scale!(y::Fun, x::Fun, α::Number)
    _check_compatible(y, x)
    cy, cx = coefficients(y), coefficients(x)
    length(cy) == length(cx) || resize!(cy, length(cx))
    cy .= scale.(cx, α)
    return y
end
function VectorInterface.scale!!(y::Fun, x::Fun, α::Number)
    if spacescompatible(y, x) &&
            (length(coefficients(y)) == length(coefficients(x)) || _resizable(coefficients(y))) &&
            VectorInterface.promote_scale(x, α) <: scalartype(y)
        return scale!(y, x, α)
    end
    return scale(x, α)
end

VectorInterface.add(y::Fun, x::Fun, α::Number, β::Number) = scale(y, β) + scale(x, α)

function VectorInterface.add!(y::Fun, x::Fun, α::Number, β::Number)
    _check_compatible(y, x)
    cy, cx = coefficients(y), coefficients(x)
    ny, nx = length(cy), length(cx)
    if ny < nx
        resize!(cy, nx)
        for i in (ny + 1):nx
            @inbounds cy[i] = scale(cx[i], α)
        end
        add!(view(cy, 1:ny), view(cx, 1:ny), α, β)
    else
        add!(view(cy, 1:nx), cx, α, β)
        β === One() || scale!(view(cy, (nx + 1):ny), β)
    end
    return y
end
function VectorInterface.add!!(y::Fun, x::Fun, α::Number, β::Number)
    if spacescompatible(y, x) &&
            (length(coefficients(y)) >= length(coefficients(x)) || _resizable(coefficients(y))) &&
            VectorInterface.promote_add(y, x, α, β) <: scalartype(y)
        return add!(y, x, α, β)
    end
    return add(y, x, α, β)
end

VectorInterface.inner(f::Fun, g::Fun) = innerproduct(f, g)

end
