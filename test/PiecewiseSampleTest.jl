using Test
using ApproxFun

@testset "piecewise sample (#635)" begin
    f = abs(Fun(sin, -5..5))
    @test integrate(f)(-4.0) ≈ -(cos(-4.0) - cos(-5.0))
    @test -(cos(-π) - cos(-5.0)) + cos(-3.0) - cos(-π) ≈ integrate(f)(-3.0)
    r = ApproxFun.sample(f,10)
    @test maximum(r) ≤ 5
    @test minimum(r) ≥ -5
end
