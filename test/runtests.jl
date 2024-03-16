versioninfo()

using ApproxFun,Base.Test

println("Domain tests")

@test !in(0.45-0.65im,Interval())
@test cumsum(ApproxFun.Flatten(([3],ApproxFun.repeated(2)))).it[2]==ApproxFun.Count(5,2)

include("MatrixTest.jl")


println("Fun tests")
include("ChebyshevTest.jl")
include("FourierTest.jl")
include("ComplexIFunTest.jl")
include("NumberTypeTest.jl")
println("Operator tests")
include("OperatorTest.jl")
println("ODE tests")
include("ODETest.jl")
println("Vector tests")
include("VectorTest.jl")
println("Singularities tests")
include("SingularitiesTest.jl")
println("Integral Equations tests")
include("IntegralEquationsTest.jl")
println("Spaces tests")
include("SpacesTest.jl")
include("JacobiTest.jl")


println("Extras tests")
include("ETDRK4Test.jl")
include("ExtrasTest.jl")
include("FractionalTest.jl")

println("Multivariate tests")
include("MultivariateTest.jl")
println("PDE tests")
include("PDETest.jl")
