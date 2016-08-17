# SplineSpace represents a Spline, right now piecewise constant HeavisideSpace is only implemented case
immutable SplineSpace{order,T} <: RealUnivariateSpace{PiecewiseInterval{T}}
    domain::PiecewiseInterval{T}
end

@compat (::Type{SplineSpace{m}}){m,T}(d::PiecewiseInterval{T}) = SplineSpace{m,T}(d)
@compat (::Type{SplineSpace{m}}){m}(d::AbstractVector) = SplineSpace{m}(PiecewiseInterval(sort(d)))

typealias HeavisideSpace{T} SplineSpace{0,T}
dimension{λ}(h::SplineSpace{λ}) = length(h.domain.points)+λ-1

Base.convert(::Type{HeavisideSpace},d::PiecewiseInterval)=HeavisideSpace{eltype(d)}(d)

Base.convert(::Type{HeavisideSpace},d::AbstractVector) =
    HeavisideSpace(PiecewiseInterval(sort(d)))

spacescompatible{λ}(a::SplineSpace{λ},b::SplineSpace{λ}) = domainscompatible(a,b)
canonicalspace(sp::HeavisideSpace) = PiecewiseSpace(map(Chebyshev,pieces(domain(sp))))


function evaluate{T<:Real}(f::Fun{HeavisideSpace{T}},x::Real)
    p = domain(f).points
    c = f.coefficients
    for k=1:length(p)-1
        if p[k] ≤ x ≤ p[k+1]
            return c[k]
        end
    end
    return zero(T)
end


function evaluate{T<:Real}(f::Fun{SplineSpace{1,T}},x::Real)
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
        function transform(S::SplineSpace{$λ},vals::Vector,plan...)
            @assert length(vals) ≤ dimension(S)
            vals
        end
        itransform(S::SplineSpace{$λ},cfs::Vector,plan...) = pad(cfs,dimension(S))
    end
end

conversion_rule{k,PS<:PolynomialSpace}(sp::HeavisideSpace,sp2::PiecewiseSpace{NTuple{k,PS}}) = sp


Conversion{kk,CC<:PolynomialSpace,DD}(a::HeavisideSpace,b::PiecewiseSpace{NTuple{kk,CC},RealBasis,DD,1}) =
    ConcreteConversion(a,b)
bandinds{HS<:HeavisideSpace,CC<:PolynomialSpace,DD,kk}(::ConcreteConversion{HS,PiecewiseSpace{NTuple{kk,CC},RealBasis,DD,1}}) =
    0,0
#bandinds{HS<:HeavisideSpace,DD,D}(::ConcreteConversion{PiecewiseSpace{ChebyshevDirichlet{1,1,D},RealBasis,DD,1},HS})=0,0
getindex{HS<:HeavisideSpace,CC<:PolynomialSpace,DD,kk}(C::ConcreteConversion{HS,PiecewiseSpace{NTuple{kk,CC},RealBasis,DD,1}},k::Integer,j::Integer) =
    k ≤ dimension(domainspace(C)) && j==k? one(eltype(C)) : zero(eltype(C))


bandinds{HS<:HeavisideSpace}(D::ConcreteDerivative{HS})=-1,0

rangespace{HS<:HeavisideSpace}(D::ConcreteDerivative{HS})=DiracSpace(domain(D).points)

function getindex{HS<:HeavisideSpace}(D::ConcreteDerivative{HS},k::Integer,j::Integer)
    n=numpieces(domain(D))
    if k≤n && j==k
        one(eltype(D))
    elseif j≤n && j==k-1
        -one(eltype(D))
    else
        zero(eltype(D))
    end
end

Base.sum{HS<:HeavisideSpace}(f::Fun{HS}) = dotu(f.coefficients,diff(space(f).domain.points))
function Base.sum{T}(f::Fun{SplineSpace{1,T}})
    vals=pad(f.coefficients,dimension(space(f)))
    dfs=diff(space(f).domain.points)
    ret=vals[1]*dfs[1]/2
    for k=2:length(vals)-1
        ret+=vals[k]*(dfs[k]+dfs[k+1])/2
    end
    ret+=vals[end]*dfs[end]/2
    ret
end


differentiate{T}(f::Fun{SplineSpace{1,T}}) =
    Fun(diff(pad(f.coefficients,dimension(space(f))))./diff(space(f).domain.points),
        HeavisideSpace(space(f).domain))

integrate{T}(f::Fun{HeavisideSpace{T}}) =
    Fun([0;cumsum(f.coefficients).*diff(space(f).domain.points)],
        SplineSpace{1,T}(space(f).domain))
