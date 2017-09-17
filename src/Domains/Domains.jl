include("Segment.jl")
include("PeriodicInterval.jl")
include("Ray.jl")
include("Circle.jl")
include("Line.jl")
include("Arc.jl")

include("UnionDomain.jl")
include("PiecewiseSegment.jl")
include("Curve.jl")

include("Point.jl")


const AffineDomain = Union{Segment,PeriodicInterval,Ray,Line}


points(d::ClosedInterval,n) = points(Domain(d),n)

# These are needed for spaces to auto-convert [a,b] to Segment
function convert(::Type{Domain},d::ClosedInterval)
    a,b=d.left,d.right
    if isinf(norm(a)) && isinf(norm(b))
        Line(d)
    elseif isinf(norm(a)) || isinf(norm(b))
        Ray(d)
    else
        Segment(d)
    end
end

# These are needed for spaces to auto-convert [a,b] to Interval
function convert(::Type{PeriodicDomain},d::ClosedInterval)
    a,b=d.left,d.right
    if isinf(norm(a)) && isinf(norm(b))
        PeriodicLine(d)
    elseif isinf(norm(a)) || isinf(norm(b))
        error("PeriodicRay not implemented")
    else
        PeriodicInterval(d)
    end
end

convert(::Type{Space},d::ClosedInterval) = Space(Domain(d))


#issubset between domains

Base.issubset(a::PeriodicInterval,b::Segment) = Segment(a.a,a.b)⊆b
Base.issubset(a::Segment,b::PeriodicInterval) = PeriodicInterval(a.a,a.b)⊆b
Base.issubset(a::Segment{T},b::PiecewiseSegment{T}) where {T<:Real} =
    a⊆Segment(first(b.points),last(b.points))
Base.issubset(a::Segment,b::Line) = first(a)∈b && last(a)∈b


function Base.intersect(a::Segment,b::Line)
    @assert a ⊆ b
    a
end

Base.intersect(b::Line,a::Segment) = intersect(a,b)


function Base.setdiff(b::Line,a::Segment)
    @assert a ⊆ b
    if first(a)>last(a)
        b\reverse(a)
    else
        Ray([first(b),first(a)]) ∪ Ray([last(a),last(b)])
    end
end

function Base.setdiff(b::Segment,a::Point)
    if !(a ⊆ b)
        b
    elseif first(b) == a.x  || last(b) == a.x
        b
    else
        Segment(first(b),a.x) ∪ Segment(a.x,last(b))
    end
end

# sort

Base.isless(d1::Segment{T1},d2::Ray{false,T2}) where {T1<:Real,T2<:Real} = d1 ≤ d2.center
Base.isless(d2::Ray{true,T2},d1::Segment{T1}) where {T1<:Real,T2<:Real} = d2.center ≤ d1


# ^
*(a::ClosedInterval,b::Domain) = Domain(a)*b
*(a::Domain,b::ClosedInterval) = a*Domain(b)

#union
Base.union(a::ClosedInterval,b::Domain) = union(Domain(a),b)
Base.union(a::Domain,b::ClosedInterval) = union(a,Domain(b))


## set minus
function Base.setdiff(d::AffineDomain,ptsin::UnionDomain{AS}) where {AS <: AbstractVector{P}} where {P <: Point}
    pts=Number.(elements(ptsin))
    isempty(pts) && return d
    tol=sqrt(eps(arclength(d)))
    da=first(d)
    isapprox(da,pts[1];atol=tol) && shift!(pts)
    isempty(pts) && return d
    db=last(d)
    isapprox(db,pts[end];atol=tol) && pop!(pts)

    sort!(pts)
    d.a > d.b && reverse!(pts)
    filter!(p->p ∈ d,pts)

    isempty(pts) && return d
    length(pts) == 1 && return d \ pts[1]

    ret = Array{Domain}(length(pts)+1)
    ret[1] = Domain(d.a..pts[1])
    for k = 2:length(pts)
        ret[k] = Domain(pts[k-1]..pts[k])
    end
    ret[end] = Domain(pts[end]..d.b)
    UnionDomain(ret)
end

function Base.setdiff(d::IntervalDomain,p::Point)
    x = Number(p)
    (x ∉ d || x ≈ first(d) || x ≈ last(d)) && return d
    DifferenceDomain(d,p)
end



# multivariate domainxs

include("multivariate.jl")
include("Disk.jl")
