

println("Fun tests")
include("IFunTest.jl")
include("FFunTest.jl")
include("ComplexIFunTest.jl")
println("ODE tests")
include("ODETest.jl")
include("VectorTest.jl")
include("OperatorTest.jl")
println("Spaces tests")
include("SingularitiesTest.jl")
include("MultivariateTest.jl")
include("SpacesTest.jl")
if isdir(Pkg.dir("FastGaussQuadrature"))
    include("JacobiTest.jl")
else
    warn("Unable to do JacobiTest.jl since FastGaussQuadrature.jl is not installed")
end
println("Speed tests")
include("SpeedTest.jl")
println("Example tests")
if isdir(Pkg.dir("Gadfly"))
    include("ExamplesTest.jl")
    include("ReadmeTest.jl")
else
    warn("Unable to do Examples since Gadfly.jl is not installed")
end