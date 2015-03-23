

immutable ContinuousSpace <: FunctionSpace{RealBasis,UnionDomain}
    domain::PiecewiseInterval
end

Space(d::PiecewiseInterval)=ContinuousSpace(d)

spacescompatible(a::ContinuousSpace,b::ContinuousSpace)=domainscompatible(a,b)
conversion_rule(a::ContinuousSpace,b::PiecewiseSpace{ChebyshevDirichlet{1,1},RealBasis})=a

function transform(S::ContinuousSpace,vals::Vector)
    n=length(vals)
    d=domain(S)
    K=length(d)
   k=div(n,K)

    PT=promote_type(eltype(d),eltype(vals))
    if k==0
        vals
    else
        ret=Array(PT,n-K+1)
        r=n-K*k
        j=1
        cfs=transform(ChebyshevDirichlet{1,1}(d[j]),vals[(j-1)*(k+1)+1:j*(k+1)])
        ret[1]=cfs[1]-cfs[2]
        ret[2]=cfs[1]+cfs[2]
        ret[K+2:K:end]=cfs[3:end]
        for j=2:r
            cfs=transform(ChebyshevDirichlet{1,1}(d[j]),vals[(j-1)*(k+1)+1:j*(k+1)])
            ret[j+1]=cfs[1]+cfs[2]
            ret[K+j+1:K:end]=cfs[3:end]
        end

        for j=r+1:K
            cfs=transform(ChebyshevDirichlet{1,1}(d[j]),vals[r*(k+1)+(j-r-1)*k+1:r*(k+1)+(j-r)*k])
            ret[j+1]=cfs[1]+cfs[2]
            ret[K+j+1:K:end]=cfs[3:end]
        end

        ret
    end
end

canonicalspace(S::ContinuousSpace)=PiecewiseSpace(map(ChebyshevDirichlet{1,1},pieces(domain(S))))




## Conversion

function addentries!{T}(C::Conversion{PiecewiseSpace{ChebyshevDirichlet{1,1},RealBasis},ContinuousSpace,T},A,kr::Range)
    d=domain(rangespace(C))
    K=length(d)
    for k=kr
        if k==1
            A[k,1]+=1
            A[k,K+1]-=1
        elseif k≤K+1
            A[k,k-1]+=1
            A[k,K+k-1]+=1
        else #K+1->
            A[k,k+K-1]+=1
        end
    end
    A
end
bandinds(C::Conversion{PiecewiseSpace{ChebyshevDirichlet{1,1},RealBasis},ContinuousSpace})=-1,length(domain(rangespace(C)))


function addentries!{T}(C::Conversion{ContinuousSpace,PiecewiseSpace{ChebyshevDirichlet{1,1},RealBasis},T},A,kr::Range)
    d=domain(rangespace(C))
    K=length(d)
    for k=kr
        if k≤K
            A[k,k]+=0.5
            A[k,k+1]+=0.5
        elseif K+1≤k≤2K
            A[k,k-K]+=-0.5
            A[k,k-K+1]+=0.5
        else #K+1->
            A[k,k-K+1]+=1
        end
    end
    A
end
bandinds(C::Conversion{ContinuousSpace,PiecewiseSpace{ChebyshevDirichlet{1,1},RealBasis}})=-length(domain(rangespace(C))),1
