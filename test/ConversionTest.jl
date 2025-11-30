using Test
using ApproxFun

@testset "Conversion" begin
    f=Fun(t->[cos(t) 0;sin(t) 1],-π..π)
    g=Fun(f,Space(PeriodicSegment(-π,π)))
    @test g(.1) ≈ f(.1)
end
