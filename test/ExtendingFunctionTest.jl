using Test
using ApproxFun

@testset "Extending function" begin
    Γ = Segment(-im,1.0-im) ∪ Curve(Fun(x->exp(0.8im)*(x+x^2-1+im*(x-4x^3+x^4)/6))) ∪ Circle(2.0,0.2)
    @test isempty(component(Γ,1)\component(Γ,1))
    @test Γ \ component(Γ,1) == component(Γ,2) ∪ component(Γ,3)
    @test norm(Fun(ones(component(Γ,1)),Γ) - Fun(x->x ∈ component(Γ,1) ? 1.0 : 0.0,Γ)) == 0
end
