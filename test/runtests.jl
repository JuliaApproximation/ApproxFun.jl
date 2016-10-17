versioninfo()

using ApproxFun,Base.Test

println("Domain tests")

@test !in(0.45-0.65im,Interval())
@test cumsum(ApproxFun.Flatten(([3],ApproxFun.repeated(2)))).it[2]==ApproxFun.Count(5,2)

@time include("MatrixTest.jl")


println("Fun tests")
@time include("ChebyshevTest.jl")
@time include("FourierTest.jl")
@time include("ComplexIFunTest.jl")
@time include("NumberTypeTest.jl")
println("Operator tests")
@time include("OperatorTest.jl")
println("ODE tests")
@time include("ODETest.jl")
println("Vector tests")
@time include("VectorTest.jl")
println("Singularities tests")
@time include("SingularitiesTest.jl")
println("Integral Equations tests")
@time include("IntegralEquationsTest.jl")
println("Spaces tests")
@time include("SpacesTest.jl")
println("Jacobi tests")
@time include("JacobiTest.jl")


println("Extras tests")
@time include("ETDRK4Test.jl")
@time include("ExtrasTest.jl")
@time include("FractionalTest.jl")

println("Multivariate tests")
@time include("MultivariateTest.jl")
println("PDE tests")
@time include("PDETest.jl")
