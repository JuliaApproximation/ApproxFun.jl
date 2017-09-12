versioninfo()

using ApproxFun, Base.Test
    import ApproxFun: Infinity, ∞

println("Helper tests")
@test ApproxFun.interlace!([-1.0],0) == [-1.0]
@test ApproxFun.interlace!([1.0,2.0],0) == [2.0,1.0]
@test ApproxFun.interlace!([1,2,3],0) == [2,1,3]
@test ApproxFun.interlace!([1,2,3,4],0) == [3,1,4,2]

@test ApproxFun.interlace!([-1.0],1) == [-1.0]
@test ApproxFun.interlace!([1.0,2.0],1) == [1.0,2.0]
@test ApproxFun.interlace!([1,2,3],1) == [1,3,2]
@test ApproxFun.interlace!([1,2,3,4],1) == [1,3,2,4]

@test exp(im*π/4)*∞ == Inf+im*Inf
@test exp(im*π/4)+∞ == ∞
@test ∞ ≠ 1

@test maximum([1,∞]) == ∞
@test minimum([1,∞]) == 1

@test (1:∞) ∩ (2:10) == (2:10)
@test (2:10) ∩ (1:∞) == (2:10)
@test (3:∞) ∩ (2:2:10) == (4:2:10)
@test (2:2:10) ∩ (3:∞) == (4:2:10)
@test (3:3:∞) ∩ (2:2:10) == 6:6:6
@test (2:2:10) ∩ (3:3:∞) == 6:6:6

@test Infinity(true)+Infinity(true) == Infinity(true)
@test Infinity(false)+Infinity(false) == Infinity(false)
@test Infinity(true)+1 == Infinity(true)
@test Infinity(false)+1 == Infinity(false)

@test ApproxFun.interlace(collect(6:10),collect(1:5)) == ApproxFun.interlace!(collect(1:10),0)
@test ApproxFun.interlace(collect(1:5),collect(6:10)) == ApproxFun.interlace!(collect(1:10),1)

@test maximum(ApproxFun.repeated(1)) == 1
@test minimum(1:∞) == 1
@test minimum(ApproxFun.flatten(([2.0],1:ApproxFun.∞))) == 1

cumsum(ApproxFun.repeated(true)) == 1:ApproxFun.∞
cumsum(ApproxFun.repeated(2)) == 2:2:ApproxFun.∞

@test 2*(1:∞) == 2:2:∞
@test 2+(1:∞) == 3:∞

# BlockInterlacer tests
@test cache(ApproxFun.BlockInterlacer((1:∞,[2],[2])))[1:6] ==
    [(1,1),(2,1),(2,2),(3,1),(3,2),(1,2)]

@test collect(ApproxFun.BlockInterlacer(([2],[2],[2]))) ==
    [(1,1),(1,2),(2,1),(2,2),(3,1),(3,2)]


# TODO: Tensorizer tests

println("Domain tests")

@test 0.45-0.65im ∉ Interval()
@test cumsum(ApproxFun.Flatten(([3],ApproxFun.repeated(2)))).it[2] ==
    ApproxFun.Count(5,2)
@test reverse(Arc(1,2,(0.1,0.2))) == Arc(1,2,(0.2,0.1))
@test 0.1 ∈ PeriodicInterval(2π,0)
@test 100.0 ∈ PeriodicInterval(0,2π)
@test -100.0 ∈ PeriodicInterval(0,2π)

@test 10.0 ∈ PeriodicLine()
@test -10.0 ∈ PeriodicLine()
@test -10.0+im ∉ PeriodicLine()

@test ApproxFun.Vec(0,0.5) ∈ PeriodicInterval(ApproxFun.Vec(0.0,0), ApproxFun.Vec(0,1))

@test ApproxFun.Vec(1,0) ∈ Circle((0.,0.),1.)

@test isambiguous(convert(ApproxFun.Point,ApproxFun.AnyDomain()))
@test isambiguous(ApproxFun.Point(ApproxFun.AnyDomain()))

@test ApproxFun.AnySegment() == ApproxFun.AnySegment()
@test ApproxFun.Point(NaN) == ApproxFun.Point(NaN)


@test ApproxFun.dimension(Domain{Float64}) == 1
@test ApproxFun.dimension(Segment{Float64}) == 1
@test ApproxFun.dimension(Interval()) == 1
@test ApproxFun.dimension(Interval()^2) == 2
@test ApproxFun.dimension(Interval()^3) == 3

@time include("MatrixTest.jl")


println("Fun tests")
@time include("ClenshawTest.jl")
@time include("ChebyshevTest.jl")
@time include("FourierTest.jl")
@time include("ComplexTest.jl")
@time include("NumberTypeTest.jl")
@time include("broadcastingtest.jl")
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

println("Multivariate tests")
@time include("MultivariateTest.jl")
println("PDE tests")
@time include("PDETest.jl")
