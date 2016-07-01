immutable PiecewiseInterval{T} <: UnivariateDomain{T}
    points::Vector{T}
    PiecewiseInterval(d::Vector{T})=new(d)
end
PiecewiseInterval(d::AbstractVector)=PiecewiseInterval{eltype(d)}(d)
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
numpieces(d::PiecewiseInterval)=length(d.points)-1
pieces(d::PiecewiseInterval,k)=d[k]
pieces{T}(d::PiecewiseInterval{T})=Interval{T}[d[k] for k=1:numpieces(d)]

for OP in (:arclength,:complexlength)
    @eval $OP(d::PiecewiseInterval)=mapreduce($OP,+,pieces(d))
end

Base.getindex(d::PiecewiseInterval,j::Integer)=Interval(d.points[j],d.points[j+1])
isperiodic(d::PiecewiseInterval)=first(d.points)==last(d.points)

Base.reverse(d::PiecewiseInterval)=PiecewiseInterval(reverse(d.points))

isambiguous(d::PiecewiseInterval)=isempty(d.points)
Base.convert{T<:Number}(::Type{PiecewiseInterval{T}},::AnyDomain)=PiecewiseInterval{T}([])
Base.convert{IT<:PiecewiseInterval}(::Type{IT},::AnyDomain)=PiecewiseInterval(Float64[])


function points(d::PiecewiseInterval,n)
   k=div(n,numpieces(d))
    r=n-numpieces(d)*k

    eltype(d)[vcat([points(d[j],k+1) for j=1:r]...);
        vcat([points(d[j],k) for j=r+1:numpieces(d)]...)]
end



Base.rand(d::PiecewiseInterval)=rand(d[rand(1:numpieces(d))])
checkpoints{T}(d::PiecewiseInterval{T})=mapreduce(checkpoints,union,pieces(d))

for OP in (:(Base.first),:(Base.last))
    @eval $OP(d::PiecewiseInterval)=$OP(d.points)
end


# Comparison with UnionDomain
for OP in (:(Base.isapprox),:(==))
    @eval begin
        $OP(a::PiecewiseInterval,b::UnionDomain)=$OP(UnionDomain(pieces(a)),b)
        $OP(b::UnionDomain,a::PiecewiseInterval)=$OP(UnionDomain(pieces(a)),b)
    end
end
