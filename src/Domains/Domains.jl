include("Segment.jl")
include("PeriodicSegment.jl")
include("Ray.jl")
include("Circle.jl")
include("Line.jl")
include("Arc.jl")

include("UnionDomain.jl")
include("PiecewiseSegment.jl")
include("Curve.jl")

include("Point.jl")


const AffineDomain = Union{Domains.AbstractInterval,Segment,PeriodicSegment,Ray,Line}


# These are needed for spaces to auto-convert [a,b] to Interval
function convert(::Type{PeriodicDomain},d::ClosedInterval)
    a,b=d.left,d.right
    a,b = float(a),float(b)
    if isinf(norm(a)) && isinf(norm(b))
        PeriodicLine(d)
    elseif isinf(norm(a)) || isinf(norm(b))
        error("PeriodicRay not implemented")
    else
        PeriodicSegment(d)
    end
end

convert(::Type{Space},d::ClosedInterval) = Space(Domain(d))

#issubset between domains

issubset(a::PeriodicSegment, b::IntervalOrSegment) = Segment(endpoints(a)...)⊆b
issubset(a::IntervalOrSegment, b::PeriodicSegment) = PeriodicSegment(endpoints(a)...)⊆b
issubset(a::IntervalOrSegment{T}, b::PiecewiseSegment{T}) where {T<:Real} =
    a⊆Segment(first(b.points),last(b.points))
issubset(a::IntervalOrSegment, b::Line) = first(a)∈b && last(a)∈b
issubset(a::Ray{angle}, b::Line{angle}) where angle = first(a) ∈ b
issubset(a::Ray{true}, b::Line{false}) = true
issubset(a::Ray{false}, b::Line{true}) = true



function intersect(a::Union{Interval,Segment,Ray},b::Line)
    @assert a ⊆ b
    a
end

function union(a::Union{Interval,Segment,Ray},b::Line)
    @assert a ⊆ b
    b
end

intersect(b::Line,a::Union{Interval,Segment,Ray}) = intersect(a,b)
union(b::Line,a::Union{Interval,Segment,Ray}) = union(a,b)


function setdiff(b::Line,a::Segment)
    @assert a ⊆ b
    if first(a)>last(a)
        b\reverseorientation(a)
    else
        Ray([first(b),first(a)]) ∪ Ray([last(a),last(b)])
    end
end

function setdiff(b::Segment,a::Point)
    if !(a ⊆ b)
        b
    elseif first(b) == a.x  || last(b) == a.x
        b
    else
        Segment(first(b),a.x) ∪ Segment(a.x,last(b))
    end
end

# sort

isless(d1::IntervalOrSegment{T1},d2::Ray{false,T2}) where {T1<:Real,T2<:Real} = d1 ≤ d2.center
isless(d2::Ray{true,T2},d1::IntervalOrSegment{T1}) where {T1<:Real,T2<:Real} = d2.center ≤ d1



## set minus
function Base.setdiff(d::AffineDomain,ptsin::UnionDomain{AS}) where {AS <: AbstractVector{P}} where {P <: Point}
    pts=Number.(elements(ptsin))
    isempty(pts) && return d
    tol=sqrt(eps(arclength(d)))
    da=leftendpoint(d)
    isapprox(da,pts[1];atol=tol) && popfirst!(pts)
    isempty(pts) && return d
    db=rightendpoint(d)
    isapprox(db,pts[end];atol=tol) && pop!(pts)

    sort!(pts)
    leftendpoint(d) > rightendpoint(d) && reverse!(pts)
    filter!(p->p ∈ d,pts)

    isempty(pts) && return d
    length(pts) == 1 && return d \ pts[1]

    ret = Array{Domain}(undef, length(pts)+1)
    ret[1] = Domain(leftendpoint(d) .. pts[1])
    for k = 2:length(pts)
        ret[k] = Domain(pts[k-1]..pts[k])
    end
    ret[end] = Domain(pts[end] .. rightendpoint(d))
    UnionDomain(ret)
end

function Base.setdiff(d::SegmentDomain,p::Point)
    x = Number(p)
    (x ∉ d || x ≈ leftendpoint(d) || x ≈ rightendpoint(d)) && return d
    DifferenceDomain(d,p)
end



# multivariate domainxs

include("multivariate.jl")
include("Disk.jl")
