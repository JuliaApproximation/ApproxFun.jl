using ApproxFun, Compat.Test

@testset "BigFloat" begin
    @testset "BigFloat constructor" begin
        single_sin = Fun(sin,Interval(0.f0,1.f0))
        double_sin = Fun(sin,Interval(0.,1.))
        big_sin = Fun(sin,Interval(parse(BigFloat,"0.0"), parse(BigFloat,"1.0")))

        @test ncoefficients(single_sin) <= ncoefficients(double_sin)
        @test ncoefficients(double_sin) <= ncoefficients(big_sin)

        @test eltype(coefficients(single_sin)) == Float32
        @test eltype(coefficients(double_sin)) == Float64
        @test eltype(coefficients(big_sin)) == BigFloat


        single_double_err = coefficients(single_sin-double_sin)[1:ncoefficients(single_sin)]
        @test norm(single_double_err) < 10eps(Float32)

        single_double_err = coefficients(double_sin-big_sin)[1:ncoefficients(double_sin)]
        @test norm(single_double_err) < 10eps(Float64)
    end

    @testset "BigFloat roots" begin
        a = Fun(Segment{BigFloat}(),BigFloat[1,2,3])
        @test norm(a.(roots(a))) == 0

        a = Fun(Taylor(Circle(BigFloat)),BigFloat[0.5,2,3])
        @test norm(a.(complexroots(a)))  ≤ eps(BigFloat)
    end

    @testset "BigFloat relative tolerance bug test Issue #557" begin
        x = Fun(BigFloat(0)..BigFloat(20_000));
        ν = BigFloat(1568)
        f = x^(ν/2-1) * exp(-x/2)
        # TODO: JacobiWeight should support general types to avoid warning here
        @test f(10.0) ≈ 10.0^(ν/2-1) * exp(-10.0/2) rtol=eps()
        ex_mathematica = parse(BigFloat, "8.6494114955713766301430915207861966687081153e778")
        @test (cumsum(f)(10.0) - ex_mathematica)/ex_mathematica ≤ eps()

        x = Fun(BigFloat(0)..BigFloat(Inf))
        ν = BigFloat(1568)
        f = x^(ν/2-1) * exp(-x/2)
        @test f(10.0) ≈ 10.0^(ν/2-1) * exp(-10.0/2) rtol=eps()
        @test_skip (cumsum(f)(10.0) - ex_mathematica)/ex_mathematica ≤ eps()
    end
end
