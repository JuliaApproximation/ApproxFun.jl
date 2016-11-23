include("Interval.jl")
include("PeriodicInterval.jl")
include("Ray.jl")
include("Circle.jl")
include("Line.jl")
include("Arc.jl")

include("UnionDomain.jl")
include("PiecewiseInterval.jl")
include("Curve.jl")

include("Point.jl")


typealias AffineDomain Union{Interval,PeriodicInterval,Ray,Line}


points(d::ClosedInterval,n) = points(Domain(d),n)

# These are needed for spaces to auto-convert [a,b] to Interval
function Base.convert(::Type{Domain},d::ClosedInterval)
    a,b=d.left,d.right
    if abs(a) == Inf && abs(b) == Inf
        Line(d)
    elseif abs(a) == Inf || abs(b) == Inf
        Ray(d)
    else
        Interval(d)
    end
end

# These are needed for spaces to auto-convert [a,b] to Interval
function Base.convert(::Type{PeriodicDomain},d::ClosedInterval)
    a,b=d.left,d.right
    if abs(a) == Inf && abs(b) == Inf
        PeriodicLine(d)
    elseif abs(a) == Inf || abs(b) == Inf
        error("PeriodicRay not implemented")
    else
        PeriodicInterval(d)
    end
end

Base.convert(::Type{Space},d::ClosedInterval) = Space(Domain(d))


#issubset between domains

Base.issubset(a::PeriodicInterval,b::Interval) = Interval(a.a,a.b)⊆b
Base.issubset(a::Interval,b::PeriodicInterval) = PeriodicInterval(a.a,a.b)⊆b
Base.issubset{T<:Real}(a::Interval{T},b::PiecewiseInterval{T}) =
    a⊆Interval(first(b.points),last(b.points))
Base.issubset(a::Interval,b::Line) = first(a)∈b && last(a)∈b


function Base.intersect(a::Interval,b::Line)
    @assert a ⊆ b
    a
end

Base.intersect(b::Line,a::Interval) = intersect(a,b)


function Base.setdiff(b::Line,a::Interval)
    @assert a ⊆ b
    if first(a)>last(a)
        b\reverse(a)
    else
        Ray([first(b),first(a)]) ∪ Ray([last(a),last(b)])
    end
end

function Base.setdiff(b::Interval,a::Point)
    if !(a ⊆ b)
        a
    elseif first(b) == a.x  || last(b) == a.x
        a
    else
        Interval(first(b),a.x) ∪ Interval(a.x,last(b))
    end
end

# sort

Base.isless{T1<:Real,T2<:Real}(d1::Interval{T1},d2::Ray{false,T2}) = d1 ≤ d2.center
Base.isless{T1<:Real,T2<:Real}(d2::Ray{true,T2},d1::Interval{T1}) = d2.center ≤ d1


#union 
Base.union(a::ClosedInterval,b::ClosedInterval) = union(Domain(a),Domain(b))
Base.union(a::ClosedInterval,b) = union(Domain(a),b)
Base.union(a,b::ClosedInterval) = union(a,Domain(b))


## set minus
\(d::ClosedInterval,x) = Domain(d) \ x
\(d::Domain,x::Number) = d \ Point(x)
function \(d::AffineDomain,x::Vector)
    isempty(x) && return d
    length(x) == 1 && return d \ x[1]

    x = sort(x)
    d.a > d.b && reverse!(x)
    filter!(p->p ∈ d,x)

    ret = Array(Domain,length(x)+1)
    ret[1] = Domain(d.a..x[1])
    for k = 2:length(x)
        ret[k] = Domain(x[k-1]..x[k])
    end
    ret[end] = Domain(x[end]..d.b)
    UnionDomain(ret)
end



# multivariate domainxs

include("multivariate.jl")
include("Disk.jl")
