

immutable ContinuousSpace <: Space{RealBasis,PiecewiseInterval,1}
    domain::PiecewiseInterval
end



Space(d::PiecewiseInterval)=ContinuousSpace(d)

isperiodic(C::ContinuousSpace)=isperiodic(domain(C))

spacescompatible(a::ContinuousSpace,b::ContinuousSpace)=domainscompatible(a,b)
conversion_rule{CD<:Tuple{Vararg{ChebyshevDirichlet{1,1}}}}(a::ContinuousSpace,b::PiecewiseSpace{CD,RealBasis})=a

function transform(S::ContinuousSpace,vals::Vector)
    n=length(vals)
    d=domain(S)
    K=numpieces(d)
    k=div(n,K)

    PT=promote_type(eltype(d),eltype(vals))
    if k==0
        vals
    elseif isperiodic(d)
        ret=Array(PT,max(K,n-K))
        r=n-K*k

        for j=1:r
            cfs=transform(ChebyshevDirichlet{1,1}(d[j]),vals[(j-1)*(k+1)+1:j*(k+1)])
            if j==1
                ret[1]=cfs[1]-cfs[2]
                ret[2]=cfs[1]+cfs[2]
            elseif j < K
                ret[j+1]=cfs[1]+cfs[2]
            end
            ret[K+j:K:end]=cfs[3:end]
        end

        for j=r+1:K
            cfs=transform(ChebyshevDirichlet{1,1}(d[j]),vals[r*(k+1)+(j-r-1)*k+1:r*(k+1)+(j-r)*k])
            if length(cfs)==1 && j <K
                ret[j+1]=cfs[1]
            elseif j==1
                ret[1]=cfs[1]-cfs[2]
                ret[2]=cfs[1]+cfs[2]
            elseif j < K
                ret[j+1]=cfs[1]+cfs[2]
            end
            ret[K+j:K:end]=cfs[3:end]
        end

        ret
    else
        ret=Array(PT,n-K+1)
        r=n-K*k

        for j=1:r
            cfs=transform(ChebyshevDirichlet{1,1}(d[j]),vals[(j-1)*(k+1)+1:j*(k+1)])
            if j==1
                ret[1]=cfs[1]-cfs[2]
            end

            ret[j+1]=cfs[1]+cfs[2]
            ret[K+j+1:K:end]=cfs[3:end]
        end

        for j=r+1:K
            cfs=transform(ChebyshevDirichlet{1,1}(d[j]),vals[r*(k+1)+(j-r-1)*k+1:r*(k+1)+(j-r)*k])
            if j==1
                ret[1]=cfs[1]-cfs[2]
            end
            ret[j+1]=cfs[1]+cfs[2]
            ret[K+j+1:K:end]=cfs[3:end]
        end

        ret
    end
end

canonicalspace(S::ContinuousSpace)=PiecewiseSpace(map(ChebyshevDirichlet{1,1},pieces(domain(S))))


## pieces

Base.vec{T}(f::Fun{ContinuousSpace,T},j::Integer)=vec(Fun(f,canonicalspace(f)),j)
Base.vec{T}(f::Fun{ContinuousSpace,T})=vec(Fun(f,canonicalspace(space(f))))
pieces{T}(f::Fun{ContinuousSpace,T})=vec(f)


function points(f::Fun{ContinuousSpace})
    n=length(f)
    d=domain(f)
    K=numpieces(d)

    m=isperiodic(d)?max(K,n+2K-1):n+K
    points(f.space,m)
end

## Conversion

coefficients(cfsin::Vector,A::ContinuousSpace,B::PiecewiseSpace)=defaultcoefficients(cfsin,A,B)
bandinds{CD<:Tuple{Vararg{ChebyshevDirichlet{1,1}}},DD}(C::Conversion{PiecewiseSpace{CD,RealBasis,DD,1},ContinuousSpace})=-1,numpieces(domain(rangespace(C)))

function addentries!{T,DD,CD<:Tuple{Vararg{ChebyshevDirichlet{1,1}}}}(C::Conversion{PiecewiseSpace{CD,RealBasis,DD,1},ContinuousSpace,T},A,kr::Range,::Colon)
    d=domain(rangespace(C))
    K=numpieces(d)
    if isperiodic(d)
        for k=kr
            if k==1
                A[k,1]+=1
                A[k,K+1]-=1
            elseif k≤K
                A[k,k-1]+=1
                A[k,K+k-1]+=1
            else #K+1->
                A[k,k+K]+=1
            end
        end
    else
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
    end
    A
end

bandinds{CD<:Tuple{Vararg{ChebyshevDirichlet{1,1}}},DD}(C::Conversion{ContinuousSpace,PiecewiseSpace{CD,RealBasis,DD,1}})=isperiodic(domainspace(C))?(1-2numpieces(domain(rangespace(C))),1):(-numpieces(domain(rangespace(C))),1)
function addentries!{T,CD<:Tuple{Vararg{ChebyshevDirichlet{1,1}}},DD}(C::Conversion{ContinuousSpace,PiecewiseSpace{CD,RealBasis,DD,1},T},A,kr::Range,::Colon)
    d=domain(domainspace(C))
    K=numpieces(d)
    if isperiodic(d)
        for k=kr
            if k<K
                A[k,k]+=0.5
                A[k,k+1]+=0.5
            elseif k==K
                A[k,k]+=0.5
                A[k,1]+=0.5
            elseif K+1≤k<2K
                A[k,k-K]+=-0.5
                A[k,k-K+1]+=0.5
            elseif k==2K
                A[k,k-K]+=-0.5
                A[k,1]+=0.5
            else #K+1->
                A[k,k-K]+=1
            end
        end
    else
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
    end
    A
end
