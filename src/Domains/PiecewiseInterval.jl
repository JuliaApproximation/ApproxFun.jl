using ApproxFun
    import ApproxFun:canonicalspace,spacescompatible,points


immutable PiecewiseInterval{T<:Number} <:Domain{T}
    points::Vector{T}
end

Base.length(d::PiecewiseInterval)=length(d.points)-1

Base.getindex(d::PiecewiseInterval,j::Integer)=Interval(d.points[j],d.points[j+1])

function points(d::PiecewiseInterval,n)
   k=div(n,length(d))
    r=n-length(d)*k

    [vcat([points(d[j],k+1) for j=1:r]...);
        vcat([points(d[j],k) for j=r+1:length(d)]...)]
end



Base.rand(d::PiecewiseInterval)=rand(d[rand(1:length(d))])
checkpoints{T}(d::PiecewiseInterval{T})=union(T[checkpoints(d[k]) for k=1:length(d)])



immutable PiecewiseIntervalSpace{S<:FunctionSpace,T} <: FunctionSpace{T,UnionDomain}
    spaces::Vector{S}
    PiecewiseSpace(::AnyDomain)=new(S[S(AnyDomain())])
    PiecewiseSpace(sp::Vector{S})=new(sp)
end


d=PiecewiseInterval([1.,2.,3.])

points(d,10)

