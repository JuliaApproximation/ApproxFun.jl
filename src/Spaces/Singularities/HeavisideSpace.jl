# SplineSpace represents a Spline, right now piecewise constant HeavisideSpace is only implemented case
immutable SplineSpace{order,T} <: RealUnivariateSpace{PiecewiseInterval{T}}
    domain::PiecewiseInterval{T}
end

typealias HeavisideSpace{T} SplineSpace{0,T}

dimension(h::HeavisideSpace)=length(h.domain.points)-1

Base.convert(::Type{HeavisideSpace},d::PiecewiseInterval)=HeavisideSpace{eltype(d)}(d)

function Base.convert(::Type{HeavisideSpace},d::AbstractVector)
    d=sort(d)
    HeavisideSpace(PiecewiseInterval(d))
end

spacescompatible(::HeavisideSpace,::HeavisideSpace) = true
canonicalspace(sp::HeavisideSpace) = PiecewiseSpace(map(Chebyshev,pieces(domain(sp))))

function transform(S::HeavisideSpace,vals::Vector,plan...)
    @assert length(vals) == dimension(S)
    vals
end

conversion_rule{k,PS<:PolynomialSpace}(sp::HeavisideSpace,sp2::PiecewiseSpace{NTuple{k,PS}})=sp


Conversion{kk,CC<:PolynomialSpace,DD}(a::HeavisideSpace,b::PiecewiseSpace{NTuple{kk,CC},RealBasis,DD,1})=ConcreteConversion(a,b)
bandinds{HS<:HeavisideSpace,CC<:PolynomialSpace,DD,kk}(::ConcreteConversion{HS,PiecewiseSpace{NTuple{kk,CC},RealBasis,DD,1}})=0,0
#bandinds{HS<:HeavisideSpace,DD,D}(::ConcreteConversion{PiecewiseSpace{ChebyshevDirichlet{1,1,D},RealBasis,DD,1},HS})=0,0
getindex{HS<:HeavisideSpace,CC<:PolynomialSpace,DD,kk}(C::ConcreteConversion{HS,PiecewiseSpace{NTuple{kk,CC},RealBasis,DD,1}},k::Integer,j::Integer) =
    k ≤ dimension(domainspace(C)) && j==k? one(eltype(C)) : zero(eltype(C))

# function addentries!{HS<:HeavisideSpace,CC<:PolynomialSpace,DD}(C::ConcreteConversion{PiecewiseSpace{CC,RealBasis,DD,1},HS},A,kr::Range,::Colon)
#    for k=kr
#         A[k,k]+=1
#     end
#     A
# end


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
