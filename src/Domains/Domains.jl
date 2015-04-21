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
