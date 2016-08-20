
include("ConstantSpace.jl")

include("IntervalSpace.jl")
include("PolynomialSpace.jl")
include("PeriodicSpace.jl")


include("Modifier/Modifier.jl")

## Union

# union_rule dictates how to create a space that both spaces can be converted to
# in this case, it means
function union_rule{S1<:Tuple{Vararg{PolynomialSpace}},
                    S2<:Tuple{Vararg{PolynomialSpace}}}(s1::PiecewiseSpace{S1},s2::PiecewiseSpace{S2})
    PiecewiseSpace(map(Space,merge(domain(s1),domain(s2)).domains))
end

function union_rule{S1<:Tuple{Vararg{PolynomialSpace}}}(s1::PiecewiseSpace{S1},s2::PolynomialSpace)
    PiecewiseSpace(map(Space,merge(domain(s1),domain(s2)).domains))
end




include("Chebyshev/Chebyshev.jl")
include("Ultraspherical/Ultraspherical.jl")
include("Fourier/Fourier.jl")
include("Singularities/Singularities.jl")
include("Jacobi/Jacobi.jl")
include("Hermite/Hermite.jl")
include("CurveSpace.jl")
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
