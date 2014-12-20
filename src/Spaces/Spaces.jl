
include("IntervalSpace.jl")
include("PeriodicSpace.jl")


include("Modifier/Modifier.jl")

include("Ultraspherical/Ultraspherical.jl")
include("Fourier/FourierSpace.jl")
include("Singularities/Singularities.jl")
include("Jacobi/Jacobi.jl")
include("Mapped/Mapped.jl")

## 2D

include("Disk/DiskSpace.jl")


typealias MappedChebyshev Union(Chebyshev,MappedSpace{Chebyshev})