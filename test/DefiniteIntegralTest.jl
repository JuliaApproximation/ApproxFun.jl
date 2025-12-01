using Test
using ApproxFun

@testset "definite integral" begin
    Σ = DefiniteIntegral()
    f1 = Fun(t->cos(cos(t)),-π..π)
    f = Fun(t->cos(cos(t)),Laurent(-π..π))
    @test sum(f1) ≈ Σ*f
end
