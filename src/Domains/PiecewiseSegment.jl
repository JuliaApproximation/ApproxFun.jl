struct PiecewiseSegment{T} <: UnivariateDomain{T}
    points::Vector{T}
    (::Type{PiecewiseSegment{T}}){T}(d::Vector{T}) = new{T}(d)
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

==(a::PiecewiseSegment,b::PiecewiseSegment) = a.points==b.points


canonicaldomain(d::PiecewiseSegment)=d
ncomponents(d::PiecewiseSegment)=length(d.points)-1
component(d::PiecewiseSegment,j::Integer) = Segment(d.points[j],d.points[j+1])
components{T}(d::PiecewiseSegment{T}) = Segment{T}[component(d,k) for k=1:ncomponents(d)]

for OP in (:arclength,:complexlength)
    @eval $OP(d::PiecewiseSegment) = mapreduce($OP,+,components(d))
end


isperiodic(d::PiecewiseSegment) = first(d.points)==last(d.points)

Base.reverse(d::PiecewiseSegment) = PiecewiseSegment(reverse(d.points))

isambiguous(d::PiecewiseSegment)=isempty(d.points)
convert{T<:Number}(::Type{PiecewiseSegment{T}},::AnyDomain)=PiecewiseSegment{T}([])
convert{IT<:PiecewiseSegment}(::Type{IT},::AnyDomain)=PiecewiseSegment(Float64[])


function points(d::PiecewiseSegment,n)
   k=div(n,ncomponents(d))
    r=n-ncomponents(d)*k

    eltype(d)[vcat([points(component(d,j),k+1) for j=1:r]...);
        vcat([points(component(d,j),k) for j=r+1:ncomponents(d)]...)]
end



Base.rand(d::PiecewiseSegment) = rand(d[rand(1:ncomponents(d))])
checkpoints{T}(d::PiecewiseSegment{T}) =
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
