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
function Domain(d::ClosedInterval)
    a,b=d.left,d.right
    if isinf(norm(a)) && isinf(norm(b))
        Line(d)
    elseif isinf(norm(a)) || isinf(norm(b))
        Ray(d)
    else
        Segment(d)
    end
end
convert(::Type{D}, d::ClosedInterval) where D<:Domain = D(d)

# These are needed for spaces to auto-convert [a,b] to Interval
function PeriodicDomain(d::ClosedInterval)
    a,b=d.left,d.right
    if isinf(norm(a)) && isinf(norm(b))
        PeriodicLine(d)
    elseif isinf(norm(a)) || isinf(norm(b))
        error("PeriodicRay not implemented")
    else
        PeriodicInterval(d)
    end
end

Space(d::ClosedInterval) = Space(Domain(d))
convert(::Type{S},d::ClosedInterval) where S<:Space =
    S(d)


#issubset between domains

issubset(a::PeriodicInterval,b::Segment) = Segment(a.a,a.b)⊆b
issubset(a::Segment,b::PeriodicInterval) = PeriodicInterval(a.a,a.b)⊆b
issubset(a::Segment{T},b::PiecewiseSegment{T}) where {T<:Real} =
    a⊆Segment(first(b.points),last(b.points))
issubset(a::Segment,b::Line) = first(a)∈b && last(a)∈b


function intersect(a::Segment,b::Line)
    @assert a ⊆ b
    a
end

intersect(b::Line,a::Segment) = intersect(a,b)


function setdiff(b::Line,a::Segment)
    @assert a ⊆ b
    if first(a)>last(a)
        b\reverse(a)
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

isless(d1::Segment{T1},d2::Ray{false,T2}) where {T1<:Real,T2<:Real} = d1 ≤ d2.center
isless(d2::Ray{true,T2},d1::Segment{T1}) where {T1<:Real,T2<:Real} = d2.center ≤ d1


# ^
*(a::ClosedInterval, b::Domain) = Domain(a)*b
*(a::Domain, b::ClosedInterval) = a*Domain(b)

#union
union(a::ClosedInterval,b::Domain) = union(Domain(a),b)
union(a::Domain,b::ClosedInterval) = union(a,Domain(b))


## set minus
\(d::Domain,x::Number) = d \ Point(x)


function setdiff(d::AffineDomain,ptsin::Vector)
    pts=copy(ptsin)
    isempty(pts) && return d
    tol=sqrt(eps(arclength(d)))
    da=first(d)
    isapprox(da,pts[1];atol=tol) && popfirst!(pts)
    isempty(pts) && return d
    db=last(d)
    isapprox(db,pts[end];atol=tol) && pop!(pts)

    sort!(pts)
    d.a > d.b && reverse!(pts)
    filter!(p->p ∈ d,pts)

    isempty(pts) && return d
    length(pts) == 1 && return d \ pts[1]

    ret = Array{Domain}(undef, length(pts)+1)
    ret[1] = Domain(d.a..pts[1])
    for k = 2:length(pts)
        ret[k] = Domain(pts[k-1]..pts[k])
    end
    ret[end] = Domain(pts[end]..d.b)
    UnionDomain(ret)
end





# multivariate domainxs

include("multivariate.jl")
include("Disk.jl")
