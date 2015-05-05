immutable HeavisideSpace{V} <: RealUnivariateSpace
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

conversion_rule(sp::HeavisideSpace,sp2::PiecewiseSpace{ChebyshevDirichlet{1,1}})=sp

bandinds{HS<:HeavisideSpace,PS<:PiecewiseSpace{ChebyshevDirichlet{1,1}}}(::Conversion{HS,PS})=0,0
bandinds{HS<:HeavisideSpace,PS<:PiecewiseSpace{ChebyshevDirichlet{1,1}}}(::Conversion{PS,HS})=0,0
function addentries!{HS<:HeavisideSpace,PS<:PiecewiseSpace{ChebyshevDirichlet{1,1}}}(C::Conversion{HS,PS},A,kr::Range)
   for k=kr
        A[k,k]+=1
    end
    A
end

function addentries!{HS<:HeavisideSpace,
                     PS<:PiecewiseSpace{ChebyshevDirichlet{1,1}}}(C::Conversion{PS,HS},A,kr::Range)
   for k=kr
        A[k,k]+=1
    end
    A
end


function bandinds{HS<:HeavisideSpace}(D::Derivative{HS})
    n=length(domain(D).points)
    csp=canonicalspace(domainspace(D))
    bandinds(StrideOperator(Derivative(csp),n-2,0,1,1))
end

function addentries!{HS<:HeavisideSpace}(D::Derivative{HS},A,kr::Range)
   n=length(domain(D).points)
    csp=canonicalspace(domainspace(D))

    if 1 in kr
        A[1,1]+=-1
        A[1,2]+=1
        A[1,3]+=-1
        A[1,4]+=-1
    end

    addentries!(StrideOperator(Derivative(csp),n-2,0,1,1),A,kr)
end