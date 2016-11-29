immutable PiecewiseSegment{T} <: UnivariateDomain{T}
    points::Vector{T}
    PiecewiseSegment(d::Vector{T})=new(d)
end
PiecewiseSegment(d::AbstractVector) = PiecewiseSegment{eltype(d)}(collect(d))
PiecewiseSegment(d...) = PiecewiseSegment([d...])

function PiecewiseSegment{IT<:Segment}(pcsin::AbstractVector{IT})
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
    PiecewiseSegment(p)
end

==(a::PiecewiseSegment,b::PiecewiseSegment)=a.points==b.points


canonicaldomain(d::PiecewiseSegment)=d
numpieces(d::PiecewiseSegment)=length(d.points)-1
pieces(d::PiecewiseSegment,k)=d[k]
pieces{T}(d::PiecewiseSegment{T}) = Segment{T}[d[k] for k=1:numpieces(d)]

for OP in (:arclength,:complexlength)
    @eval $OP(d::PiecewiseSegment)=mapreduce($OP,+,pieces(d))
end

Base.getindex(d::PiecewiseSegment,j::Integer) = Segment(d.points[j],d.points[j+1])
isperiodic(d::PiecewiseSegment) = first(d.points)==last(d.points)

Base.reverse(d::PiecewiseSegment) = PiecewiseSegment(reverse(d.points))

isambiguous(d::PiecewiseSegment)=isempty(d.points)
Base.convert{T<:Number}(::Type{PiecewiseSegment{T}},::AnyDomain)=PiecewiseSegment{T}([])
Base.convert{IT<:PiecewiseSegment}(::Type{IT},::AnyDomain)=PiecewiseSegment(Float64[])


function points(d::PiecewiseSegment,n)
   k=div(n,numpieces(d))
    r=n-numpieces(d)*k

    eltype(d)[vcat([points(d[j],k+1) for j=1:r]...);
        vcat([points(d[j],k) for j=r+1:numpieces(d)]...)]
end



Base.rand(d::PiecewiseSegment) = rand(d[rand(1:numpieces(d))])
checkpoints{T}(d::PiecewiseSegment{T}) = mapreduce(checkpoints,union,pieces(d))

for OP in (:(Base.first),:(Base.last))
    @eval $OP(d::PiecewiseSegment) = $OP(d.points)
end


# Comparison with UnionDomain
for OP in (:(Base.isapprox),:(==))
    @eval begin
        $OP(a::PiecewiseSegment,b::UnionDomain) = $OP(UnionDomain(pieces(a)),b)
        $OP(b::UnionDomain,a::PiecewiseSegment) = $OP(UnionDomain(pieces(a)),b)
    end
end
