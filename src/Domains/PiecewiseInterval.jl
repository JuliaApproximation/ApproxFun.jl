using ApproxFun
    import ApproxFun:canonicalspace,spacescompatible,points,RealBasis,transform,checkpoints,pieces,canonicalspace,spacescompatible,
                    domainscompatible,conversion_rule


immutable PiecewiseInterval{T<:Number} <:Domain{T}
    points::Vector{T}
end

PiecewiseInterval(d::Number...)=PiecewiseInterval([d...])
==(a::PiecewiseInterval,b::PiecewiseInterval)=a.points==b.points

Base.length(d::PiecewiseInterval)=length(d.points)-1
Base.getindex(d::PiecewiseInterval,j::Integer)=Interval(d.points[j],d.points[j+1])

function points(d::PiecewiseInterval,n)
   k=div(n,length(d))
    r=n-length(d)*k

    [vcat([points(d[j],k+1) for j=1:r]...);
        vcat([points(d[j],k) for j=r+1:length(d)]...)]
end

pieces(d::PiecewiseInterval,k)=d[k]
pieces{T}(d::PiecewiseInterval{T})=Interval{T}[d[k] for k=1:length(d)]

Base.rand(d::PiecewiseInterval)=rand(d[rand(1:length(d))])
checkpoints{T}(d::PiecewiseInterval{T})=mapreduce(checkpoints,union,pieces(d))
Base.first(d::PiecewiseInterval)=first(d.points)


immutable ContinuousSpace <: FunctionSpace{RealBasis,UnionDomain}
    domain::PiecewiseInterval
end

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


