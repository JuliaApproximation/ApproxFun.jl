
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
#include("Mapped/Mapped.jl")



#typealias MappedChebyshev Union{Chebyshev{Interval{Float64}},MappedSpace{Chebyshev{Interval{Float64}}}}






## Derivative

function invfromcanonicalD(S::Laurent{PeriodicLine{false}})
    d=domain(S)
    @assert d.center==0  && d.L==1.0
    a=Fun([1.,.5,.5],Laurent())
end

function invfromcanonicalD(S::LaurentDirichlet{PeriodicLine{false}})
    d=domain(S)
    @assert d.center==0  && d.L==1.0
    a=Fun([1.,.5,.5],Laurent())
end
