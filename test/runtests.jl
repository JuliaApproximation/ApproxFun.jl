versioninfo()

using ApproxFun,Base.Test


println("Helper tests")
@test ApproxFun.interlace!([-1.0],1) == [-1.0]
@test ApproxFun.interlace!([1.0,2.0],1) == [1.0,2.0]
@test ApproxFun.interlace!([1,2,3],1) == [1,3,2]
@test ApproxFun.interlace!([1,2,3,4],1) == [1,3,2,4]

import ApproxFun: Infinity, ∞

@test exp(im*π/4)*∞ == Inf+im*Inf
@test exp(im*π/4)+∞ == ∞
@test ∞ ≠ 1

@test maximum([1,∞]) == ∞
@test minimum([1,∞]) == 1

@test Infinity(true)+Infinity(true) == Infinity(true)
@test Infinity(false)+Infinity(false) == Infinity(false)
@test Infinity(true)+1 == Infinity(true)
@test Infinity(false)+1 == Infinity(false)

@test ApproxFun.interlace(collect(6:10),collect(1:5)) == ApproxFun.interlace!(collect(1:10),0)
@test ApproxFun.interlace(collect(1:5),collect(6:10)) == ApproxFun.interlace!(collect(1:10),1)


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
