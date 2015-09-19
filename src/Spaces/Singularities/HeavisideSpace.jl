immutable HeavisideSpace{V} <: RealUnivariateSpace{AnyDomain}
    domain::V
    HeavisideSpace(d::V)=new(d)
end

HeavisideSpace(d::PiecewiseInterval)=HeavisideSpace{typeof(d)}(d)

function HeavisideSpace(d::Vector)
    d=sort(d)
   HeavisideSpace(PiecewiseInterval(d))
end

spacescompatible(::HeavisideSpace,::HeavisideSpace)=true
canonicalspace(sp::HeavisideSpace)=PiecewiseSpace(map(ChebyshevDirichlet{1,1},pieces(domain(sp))))

conversion_rule{D}(sp::HeavisideSpace,sp2::PiecewiseSpace{ChebyshevDirichlet{1,1,D}})=sp

bandinds{HS<:HeavisideSpace,D}(::Conversion{HS,PiecewiseSpace{ChebyshevDirichlet{1,1,D},RealBasis,1}})=0,0
bandinds{HS<:HeavisideSpace,D}(::Conversion{PiecewiseSpace{ChebyshevDirichlet{1,1,D},RealBasis,1},HS})=0,0
function addentries!{HS<:HeavisideSpace,D}(C::Conversion{HS,PiecewiseSpace{ChebyshevDirichlet{1,1,D},RealBasis,1}},A,kr::Range)
   for k=kr
        A[k,k]+=1
    end
    A
end

function addentries!{HS<:HeavisideSpace,D}(C::Conversion{PiecewiseSpace{ChebyshevDirichlet{1,1,D},RealBasis,1},HS},A,kr::Range)
   for k=kr
        A[k,k]+=1
    end
    A
end


function bandinds{HS<:HeavisideSpace}(D::Derivative{HS})
    n=length(domain(D).points)
    csp=canonicalspace(domainspace(D))
    bi=bandinds(StrideOperator(Derivative(csp),n-2,0,1,1))
    bi[1],max(bi[2],n)
end

function rangespace{HS<:HeavisideSpace}(D::Derivative{HS})
    csp=canonicalspace(domainspace(D))
    d=domain(D)
    DiracSpace(rangespace(Derivative(csp)),d.points[2:end-1])
end

function addentries!{HS<:HeavisideSpace}(D::Derivative{HS},A,kr::Range)
   n=length(domain(D).points)
    csp=canonicalspace(domainspace(D))

    for k=kr∩(1:(n-2))
        A[k,k]+=-1
        A[k,k+1]+= 1
        A[k,k+n-1]+=-1
        A[k,k+n]+=-1
    end


    addentries!(StrideOperator(Derivative(csp),n-2,0,1,1),A,kr)
end
