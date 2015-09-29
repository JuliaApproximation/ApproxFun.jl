
include("ConstantSpace.jl")

include("IntervalSpace.jl")
include("PolynomialSpace.jl")
include("PeriodicSpace.jl")


include("Modifier/Modifier.jl")

include("Ultraspherical/Ultraspherical.jl")
include("Fourier/Fourier.jl")
include("Singularities/Singularities.jl")
include("Jacobi/Jacobi.jl")
include("Hermite/Hermite.jl")
include("Mapped/Mapped.jl")



typealias MappedChebyshev Union{Chebyshev{Interval{Float64}},MappedSpace{Chebyshev{Interval{Float64}}}}
