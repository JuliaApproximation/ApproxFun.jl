
include("IntervalDomainSpace.jl")
include("PeriodicDomainSpace.jl")


include("Modifier/Modifier.jl")

include("Ultraspherical/UltrasphericalSpace.jl")
include("Fourier/FourierSpace.jl")
include("Singularities/Singularities.jl")
include("Jacobi/Jacobi.jl")
include("Mapped/Mapped.jl")

## 2D

include("Disk/DiskSpace.jl")


typealias MappedChebyshevSpace Union(ChebyshevSpace,MappedSpace{ChebyshevSpace})