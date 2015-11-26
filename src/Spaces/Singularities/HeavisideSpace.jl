# SplineSpace represents a Spline, right now piecewise constant HeavisideSpace is only implemented case
immutable SplineSpace{order,T} <: RealUnivariateSpace{PiecewiseInterval{T}}
    domain::PiecewiseInterval{T}
end

typealias HeavisideSpace{T} SplineSpace{0,T}

dimension(h::HeavisideSpace)=length(h.domain.points)-1

Base.convert(::Type{HeavisideSpace},d::PiecewiseInterval)=HeavisideSpace{eltype(d)}(d)

function Base.convert(::Type{HeavisideSpace},d::Vector)
    d=sort(d)
   HeavisideSpace(PiecewiseInterval(d))
end

spacescompatible(::HeavisideSpace,::HeavisideSpace)=true
canonicalspace(sp::HeavisideSpace)=PiecewiseSpace(map(Chebyshev,pieces(domain(sp))))


conversion_rule{k,PS<:PolynomialSpace}(sp::HeavisideSpace,sp2::PiecewiseSpace{NTuple{k,PS}})=sp

bandinds{HS<:HeavisideSpace,CC<:PolynomialSpace,DD,kk}(::Conversion{HS,PiecewiseSpace{NTuple{kk,CC},RealBasis,DD,1}})=0,0
#bandinds{HS<:HeavisideSpace,DD,D}(::Conversion{PiecewiseSpace{ChebyshevDirichlet{1,1,D},RealBasis,DD,1},HS})=0,0
function addentries!{HS<:HeavisideSpace,CC<:PolynomialSpace,DD,kk}(C::Conversion{HS,PiecewiseSpace{NTuple{kk,CC},RealBasis,DD,1}},A,kr::Range,::Colon)
    d=dimension(domainspace(C))
    for k=kr
        k ≤ d && (A[k,k]+=1)
    end
    A
end

# function addentries!{HS<:HeavisideSpace,CC<:PolynomialSpace,DD}(C::Conversion{PiecewiseSpace{CC,RealBasis,DD,1},HS},A,kr::Range,::Colon)
#    for k=kr
#         A[k,k]+=1
#     end
#     A
# end


bandinds{HS<:HeavisideSpace}(D::Derivative{HS})=-1,0

rangespace{HS<:HeavisideSpace}(D::Derivative{HS})=DiracSpace(domain(D).points)

function addentries!{HS<:HeavisideSpace}(D::Derivative{HS},A,kr::Range,::Colon)
    n=numpieces(domain(D))

    1 in kr && (A[1,1]+=1)

    for k=kr∩(2:n)
        A[k,k-1]+=-1
        A[k,k]+= 1
    end

    n in kr && (A[n+1,n]+=-1)

    A
end
