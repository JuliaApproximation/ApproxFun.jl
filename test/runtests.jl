

println("Fun tests")
include("ChebyshevTest.jl")
include("FourierTest.jl")
include("ComplexIFunTest.jl")
include("NumberTypeTest.jl")
println("ODE tests")
include("ODETest.jl")
println("Vector tests")
include("VectorTest.jl")
println("OperatorTest")
include("OperatorTest.jl")
println("Spaces tests")
include("SingularitiesTest.jl")
include("SigmaTest.jl")
include("SpacesTest.jl")
if isdir(Pkg.dir("FastGaussQuadrature"))
    include("JacobiTest.jl")
else
    warn("Unable to do JacobiTest.jl since FastGaussQuadrature.jl is not installed")
end

println("Multivariate tests")
include("MultivariateTest.jl")
println("Speed tests")
include("SpeedTest.jl")


println("Example tests")
if isdir(Pkg.dir("Gadfly"))
    include("ExamplesTest.jl")
    include("ReadmeTest.jl")
else
    warn("Unable to do Examples since Gadfly.jl is not installed")
end