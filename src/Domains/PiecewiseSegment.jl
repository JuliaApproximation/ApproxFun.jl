struct PiecewiseSegment{T} <: UnivariateDomain{T}
    points::Vector{T}
    PiecewiseSegment{T}(d::Vector{T}) where {T} = new{T}(d)
end
PiecewiseSegment(d::AbstractVector) = PiecewiseSegment{eltype(d)}(collect(d))
PiecewiseSegment(d...) = PiecewiseSegment(collect(d))

function PiecewiseSegment(pcsin::AbstractVector{IT}) where IT<:Segment
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

==(a::PiecewiseSegment,b::PiecewiseSegment) = a.points==b.points


canonicaldomain(d::PiecewiseSegment)=d
ncomponents(d::PiecewiseSegment)=length(d.points)-1
component(d::PiecewiseSegment,j::Integer) = Segment(d.points[j],d.points[j+1])
components(d::PiecewiseSegment{T}) where {T} = Segment{T}[component(d,k) for k=1:ncomponents(d)]

for OP in (:arclength,:complexlength)
    @eval $OP(d::PiecewiseSegment) = mapreduce($OP,+,components(d))
end


isperiodic(d::PiecewiseSegment) = first(d.points)==last(d.points)

Base.reverse(d::PiecewiseSegment) = PiecewiseSegment(reverse(d.points))

isambiguous(d::PiecewiseSegment)=isempty(d.points)
convert(::Type{PiecewiseSegment{T}},::AnyDomain) where {T<:Number}=PiecewiseSegment{T}([])
convert(::Type{IT},::AnyDomain) where {IT<:PiecewiseSegment}=PiecewiseSegment(Float64[])


function points(d::PiecewiseSegment,n)
   k=div(n,ncomponents(d))
    r=n-ncomponents(d)*k

    eltype(d)[vcat([points(component(d,j),k+1) for j=1:r]...);
        vcat([points(component(d,j),k) for j=r+1:ncomponents(d)]...)]
end



Base.rand(d::PiecewiseSegment) = rand(d[rand(1:ncomponents(d))])
checkpoints(d::PiecewiseSegment{T}) where {T} =
    mapreduce(checkpoints,union,components(d))::Vector{T}

for OP in (:(Base.first),:(Base.last))
    @eval $OP(d::PiecewiseSegment) = $OP(d.points)
end


# Comparison with UnionDomain
for OP in (:(Base.isapprox),:(==))
    @eval begin
        $OP(a::PiecewiseSegment,b::UnionDomain) = $OP(UnionDomain(components(a)),b)
        $OP(b::UnionDomain,a::PiecewiseSegment) = $OP(UnionDomain(components(a)),b)
    end
end
