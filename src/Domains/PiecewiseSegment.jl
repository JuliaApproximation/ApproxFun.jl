struct PiecewiseSegment{T} <: Domain{T}
    points::Vector{T}
    PiecewiseSegment{T}(d::Vector{T}) where {T} = new{T}(d)
end
PiecewiseSegment(d::AbstractVector) = PiecewiseSegment{eltype(d)}(collect(d))
PiecewiseSegment(d...) = PiecewiseSegment(collect(mapreduce(eltype,promote_type,d),d))

function PiecewiseSegment(pcsin::AbstractVector{IT}) where IT<:IntervalOrSegment
    pcs=collect(pcsin)
    p=∂(pop!(pcs))
    successful=true
    while successful
        successful=false
        for k=1:length(pcs)
            if leftendpoint(pcs[k]) == last(p)
                push!(p,rightendpoint(pcs[k]))
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

indomain(x, d::PiecewiseSegment) = any(in.(x, components(d)))


canonicaldomain(d::PiecewiseSegment)=d
ncomponents(d::PiecewiseSegment)=length(d.points)-1
component(d::PiecewiseSegment,j::Integer) = Segment(d.points[j],d.points[j+1])
components(d::PiecewiseSegment{T}) where {T} = Segment{T}[component(d,k) for k=1:ncomponents(d)]

for OP in (:arclength,:complexlength)
    @eval $OP(d::PiecewiseSegment) = mapreduce($OP,+,components(d))
end


isperiodic(d::PiecewiseSegment) = first(d.points) == last(d.points)

reverseorientation(d::PiecewiseSegment) = PiecewiseSegment(reverse(d.points))

isambiguous(d::PiecewiseSegment) = isempty(d.points)
convert(::Type{PiecewiseSegment{T}},::AnyDomain) where {T<:Number} = PiecewiseSegment{T}([])
convert(::Type{IT},::AnyDomain) where {IT<:PiecewiseSegment}=PiecewiseSegment(Float64[])


function points(d::PiecewiseSegment,n)
   k=div(n,ncomponents(d))
    r=n-ncomponents(d)*k

    float(eltype(d))[vcat([points(component(d,j),k+1) for j=1:r]...);
        vcat([points(component(d,j),k) for j=r+1:ncomponents(d)]...)]
end



rand(d::PiecewiseSegment) = rand(d[rand(1:ncomponents(d))])
checkpoints(d::PiecewiseSegment{T}) where {T} = mapreduce(checkpoints,union,components(d))

leftendpoint(d::PiecewiseSegment) = first(d.points)
rightendpoint(d::PiecewiseSegment) = last(d.points)


# Comparison with UnionDomain
for OP in (:(isapprox),:(==))
    @eval begin
        $OP(a::PiecewiseSegment,b::UnionDomain) = $OP(UnionDomain(components(a)),b)
        $OP(b::UnionDomain,a::PiecewiseSegment) = $OP(UnionDomain(components(a)),b)
    end
end


function union(S::PiecewiseSegment{<:Real}, D::IntervalOrSegment{<:Real})
    isempty(D) && return S
    a,b = endpoints(D)
    (a ∈ S || b ∈ S) && return PiecewiseSegment(sort!(union(S.points, a, b)))
    UnionDomain(S, D)
end
union(D::IntervalOrSegment{<:Real}, S::PiecewiseSegment{<:Real}) = union(S,D)
