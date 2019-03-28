using ApproxFun, LinearAlgebra, Test
    import ApproxFun: Infinity, ∞

@testset "Helper" begin
    @testset "interlace" begin
        @test ApproxFun.interlace!([-1.0],0) == [-1.0]
        @test ApproxFun.interlace!([1.0,2.0],0) == [2.0,1.0]
        @test ApproxFun.interlace!([1,2,3],0) == [2,1,3]
        @test ApproxFun.interlace!([1,2,3,4],0) == [3,1,4,2]

        @test ApproxFun.interlace!([-1.0],1) == [-1.0]
        @test ApproxFun.interlace!([1.0,2.0],1) == [1.0,2.0]
        @test ApproxFun.interlace!([1,2,3],1) == [1,3,2]
        @test ApproxFun.interlace!([1,2,3,4],1) == [1,3,2,4]

        @test ApproxFun.interlace(collect(6:10),collect(1:5)) == ApproxFun.interlace!(collect(1:10),0)
        @test ApproxFun.interlace(collect(1:5),collect(6:10)) == ApproxFun.interlace!(collect(1:10),1)
    end

    @testset "Iterators" begin
        @test cache(ApproxFun.BlockInterlacer((1:∞,[2],[2])))[1:6] ==
            [(1,1),(2,1),(2,2),(3,1),(3,2),(1,2)]

        @test collect(ApproxFun.BlockInterlacer(([2],[2],[2]))) ==
            [(1,1),(1,2),(2,1),(2,2),(3,1),(3,2)]
    end

    # TODO: Tensorizer tests
end

@testset "Domain" begin
    @test 0.45-0.65im ∉ Segment(-1,1)

    @test reverseorientation(Arc(1,2,(0.1,0.2))) == Arc(1,2,(0.2,0.1))
    @test 0.1 ∈ PeriodicSegment(2π,0)
    @test 100.0 ∈ PeriodicSegment(0,2π)
    @test -100.0 ∈ PeriodicSegment(0,2π)


    @test ApproxFun.AnySegment() == ApproxFun.AnySegment()

    @test 10.0 ∈ PeriodicLine()
    @test -10.0 ∈ PeriodicLine()
    @test -10.0+im ∉ PeriodicLine()

    @test ApproxFun.Vec(0,0.5) ∈ PeriodicSegment(ApproxFun.Vec(0.0,0), ApproxFun.Vec(0,1))

    @test ApproxFun.dimension(Domain{Float64}) == 1
    @test ApproxFun.dimension(Segment{Float64}) == 1
    @test ApproxFun.dimension(ChebyshevInterval()) == 1
    @test ApproxFun.dimension(ChebyshevInterval()^2) == 2
    @test ApproxFun.dimension(ChebyshevInterval()^3) == 3

    @test ApproxFun.Vec(1,0) ∈ Circle((0.,0.),1.)

    @test isambiguous(convert(ApproxFun.Point,ApproxFun.AnyDomain()))
    @test isambiguous(ApproxFun.Point(ApproxFun.AnyDomain()))

    @test_skip ApproxFun.Point(NaN) == ApproxFun.Point(NaN)
end

@time include("MatrixTest.jl")
@time include("ClenshawTest.jl")
@time include("ChebyshevTest.jl")
@time include("FourierTest.jl")
@time include("ComplexTest.jl")
@time include("NumberTypeTest.jl")
@time include("broadcastingtest.jl")
@time include("OperatorTest.jl")
@time include("ODETest.jl")
@time include("EigTest.jl")
@time include("VectorTest.jl")
@time include("IntegralEquationsTest.jl")
@time include("JacobiTest.jl")
@time include("LaguerreTest.jl")
@time include("HermiteTest.jl")
@time include("ETDRK4Test.jl")
@time include("MultivariateTest.jl")
@time include("PDETest.jl")
@time include("ExtrasTest.jl")
