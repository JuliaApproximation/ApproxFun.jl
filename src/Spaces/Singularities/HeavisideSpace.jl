# SplineSpace represents a Spline, right now piecewise constant HeavisideSpace is only implemented case
struct SplineSpace{order,T,R} <: Space{PiecewiseSegment{T},R}
    domain::PiecewiseSegment{T}
end

SplineSpace{m,T}(d::PiecewiseSegment{T}) where {m,T} = SplineSpace{m,T,real(eltype(T))}(d)
SplineSpace{m,T}(d::AbstractVector) where {m,T} = SplineSpace{m}(PiecewiseSegment(sort(d)))

SplineSpace{m}(d::PiecewiseSegment{T}) where {m,T} = SplineSpace{m,T,real(eltype(T))}(d)
SplineSpace{m}(d::AbstractVector) where {m} = SplineSpace{m}(PiecewiseSegment(sort(d)))

const HeavisideSpace{T,R} = SplineSpace{0,T,R}
dimension(h::SplineSpace{λ}) where {λ} = length(h.domain.points)+λ-1

convert(::Type{HeavisideSpace},d::PiecewiseSegment) = HeavisideSpace{eltype(d)}(d)

convert(::Type{HeavisideSpace},d::AbstractVector) =
    HeavisideSpace(PiecewiseSegment(sort(d)))

spacescompatible(a::SplineSpace{λ},b::SplineSpace{λ}) where {λ} = domainscompatible(a,b)
canonicalspace(sp::HeavisideSpace) = PiecewiseSpace(map(Chebyshev,components(domain(sp))))


function evaluate(f::Fun{HeavisideSpace{T,R}},x::Real) where {T<:Real,R}
    p = domain(f).points
    c = f.coefficients
    for k=1:length(p)-1
        if p[k] ≤ x ≤ p[k+1]
            return c[k]
        end
    end
    return zero(T)
end


function evaluate(f::Fun{SplineSpace{1,T,R}},x::Real) where {T<:Real,R}
    p = domain(f).points
    c = f.coefficients
    for k=1:length(p)-1
        if p[k] ≤ x ≤ p[k+1]
            return (x-p[k])*c[k+1]/(p[k+1]-p[k]) + (p[k+1]-x)*c[k]/(p[k+1]-p[k])
        end
    end
    return zero(T)
end


function points(sp::HeavisideSpace,n)
    x=sp.domain.points
    (x[1:end-1] + diff(x)/2)[1:n]
end

points(sp::SplineSpace{1},n) = sp.domain.points[1:n]

for λ = [0,1]
    @eval begin
        function transform(S::SplineSpace{$λ},vals::AbstractVector,plan...)
            @assert length(vals) ≤ dimension(S)
            vals
        end
        itransform(S::SplineSpace{$λ},cfs::AbstractVector,plan...) = pad(cfs,dimension(S))
    end
end

conversion_rule(sp::HeavisideSpace,sp2::PiecewiseSpace{NTuple{k,PS}}) where {k,PS<:PolynomialSpace} = sp


Conversion(a::HeavisideSpace,b::PiecewiseSpace{NTuple{kk,CC},DD,RR}) where {kk,CC<:PolynomialSpace,DD<:UnivariateDomain,RR<:Real} =
    ConcreteConversion(a,b)
bandinds(::ConcreteConversion{HS,PiecewiseSpace{NTuple{kk,CC},DD,RR}}) where {HS<:HeavisideSpace,CC<:PolynomialSpace,DD<:UnivariateDomain,RR<:Real,kk} =
    0,0

getindex(C::ConcreteConversion{HS,PiecewiseSpace{NTuple{kk,CC},DD,RR}},k::Integer,j::Integer) where {HS<:HeavisideSpace,CC<:PolynomialSpace,DD<:UnivariateDomain,RR<:Real,kk} =
    k ≤ dimension(domainspace(C)) && j==k? one(eltype(C)) : zero(eltype(C))

Base.sum(f::Fun{HS}) where {HS<:HeavisideSpace} = dotu(f.coefficients,diff(space(f).domain.points))
function Base.sum(f::Fun{SplineSpace{1,T,R}}) where {T,R}
    vals=pad(f.coefficients,dimension(space(f)))
    dfs=diff(space(f).domain.points)
    ret=vals[1]*dfs[1]/2
    for k=2:length(vals)-1
        ret+=vals[k]*(dfs[k]+dfs[k+1])/2
    end
    ret+=vals[end]*dfs[end]/2
    ret
end

#diffentiate HeavisideSpace
function differentiate(f::Fun{<:HeavisideSpace})
    dp=domain(f).points
    cfs=f.coefficients
    diff=0.0
    for n=1:length(cfs)
        diff=diff+cfs[n]*(DiracDelta(dp[n])-DiracDelta(dp[n+1]))
    end
    return diff
end

#Derivative Operator for HeavisideSpace
function Derivative(H::HeavisideSpace)
    ConcreteDerivative(H)
end

bandinds(D::ConcreteDerivative{H}) where {H<:HeavisideSpace} = (-1,0)
rangespace(D::ConcreteDerivative{H}) where {H<:HeavisideSpace} = domain(D).points[end]==Inf?DiracSpace(domain(D).points[1:end-1]):DiracSpace(domain(D).points)

function getindex(D::ConcreteDerivative{H,T},k::Integer,j::Integer) where {H<:HeavisideSpace,T}
    if k==j
        one(T)
    elseif k==j+1
        -one(T)
    else
        zero(T)
    end
end
    
differentiate(f::Fun{SplineSpace{1,T,R}}) where {T,R} =
    Fun(HeavisideSpace(space(f).domain),
        diff(pad(f.coefficients,dimension(space(f))))./diff(space(f).domain.points))

integrate(f::Fun{HeavisideSpace{T,R}}) where {T,R} =
    Fun(SplineSpace{1,T,R}(space(f).domain),
        [0;cumsum(f.coefficients).*diff(space(f).domain.points)])
