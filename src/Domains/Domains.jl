include("Interval.jl")
include("PeriodicInterval.jl")
include("Ray.jl")
include("Circle.jl")
include("Line.jl")
include("Arc.jl")

include("UnionDomain.jl")
include("PiecewiseInterval.jl")
include("Curve.jl")

include("multivariate.jl")

typealias AffineDomain Union(Interval,PeriodicInterval,Ray,Line)



#TODO: in 0.4 change to Domain(d::Vector{T})
# These are needed for spaces to auto-convert [a,b] to Interval
function Base.convert{D<:Domain,T<:Number}(::Type{D},d::Vector{T})
    @assert length(d) >1

    if length(d) == 2
        if abs(d[1]) == Inf && abs(d[2]) == Inf
            Line(d)
        elseif abs(d[2]) == Inf || abs(d[1]) == Inf
            Ray(d)
        else
            Interval(d[1],d[2])
        end
    else
        # TODO: convert to PiecewiseInterval
        convert(D,d[1:2])∪convert(D,d[2:end])
    end
end

function Base.convert{D<:PeriodicDomain,T<:Number}(::Type{D},d::Vector{T})
    @assert length(d) == 2

    if abs(d[1]) ==Inf
        PeriodicLine(d)
    else
        PeriodicInterval(d[1],d[2])
    end
end


Base.promote_rule{D<:Domain,T<:Number}(::Type{D},::Type{Vector{T}})=UnivariateDomain{T}
Base.promote_rule{D<:PeriodicDomain,T<:Number}(::Type{D},::Type{Vector{T}})=PeriodicDomain{T}
