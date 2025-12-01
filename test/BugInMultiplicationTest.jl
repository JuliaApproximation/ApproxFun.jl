using Test
using ApproxFun
using ApproxFunBase

@testset "Bug in Multiplication" begin
    dom = Interval(0.001, 1) × PeriodicSegment(-pi, pi)
    @test blocklengths(Space(dom)) == 2:2:∞
    r,r2 = Fun((r,t) -> [r;r^2], dom)
    @test r(0.1,0.2) ≈ 0.1
    @test r2(0.1,0.2) ≈ 0.1^2
    sp = Space(dom)
    Dr = Derivative(sp, [1,0])
    @test blockbandwidths(Dr) == (-1,1)
    @test subblockbandwidths(Dr)  == (1,3)
    Dθ = Derivative(sp, [0,1])
    Mr = Multiplication(Fun( (r, θ) -> r, sp ), sp)
    rDr = Mr * Dr
    testbandedblockbandedoperator(rDr)
end
