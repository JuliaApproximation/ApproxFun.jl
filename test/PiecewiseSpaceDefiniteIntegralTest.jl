using Test
using ApproxFun
using Random

@testset "Piecewise space definite integral" begin
    Γ=Segment(-im,1.0-im) ∪ Curve(Fun(x->exp(0.8im)*(x+x^2-1+im*(x-4x^3+x^4)/6))) ∪ Circle(2.0,0.2)
    z=Fun(Γ)
    S=PiecewiseSpace(map(d->isa(d,Circle) ? Fourier(d) : JacobiWeight(0.5,0.5,Ultraspherical(1,d)),components(Γ)))
    B=DefiniteLineIntegral(S)
    Random.seed!(0)
    f=Fun(S,rand(20))
    @test B*f ≈ linesum(component(f,1)) + linesum(component(f,2)) + linesum(component(f,3))
end
