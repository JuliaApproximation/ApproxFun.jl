immutable PiecewiseInterval{T<:Number} <: UnivariateDomain{T}
    points::Vector{T}
end

PiecewiseInterval(d::Number...)=PiecewiseInterval([d...])

function PiecewiseInterval{IT<:Interval}(pcsin::AbstractVector{IT})
    pcs=collect(pcsin)
    p=∂(pop!(pcs))
    successful=true
    while successful
        successful=false
        for k=1:length(pcs)
            if first(pcs[k])==last(p)
                push!(p,last(pcs[k]))
                deleteat!(pcs,k)
                successful=true
                break
            end
        end
    end
    @assert isempty(pcs)
    PiecewiseInterval(p)
end

==(a::PiecewiseInterval,b::PiecewiseInterval)=a.points==b.points


canonicaldomain(d::PiecewiseInterval)=d

Base.length(d::PiecewiseInterval)=length(d.points)-1
Base.getindex(d::PiecewiseInterval,j::Integer)=Interval(d.points[j],d.points[j+1])
isperiodic(d::PiecewiseInterval)=first(d.points)==last(d.points)

isambiguous(d::PiecewiseInterval)=isempty(d.points)
Base.convert{T<:Number}(::Type{PiecewiseInterval{T}},::AnyDomain)=PiecewiseInterval{T}([])
Base.convert{IT<:PiecewiseInterval}(::Type{IT},::AnyDomain)=PiecewiseInterval(Float64[])


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

for OP in (:(Base.first),:(Base.last))
    @eval $OP(d::PiecewiseInterval)=$OP(d.points)
end
