using Test
using ApproxFun

@testset "Null space" begin
    d=ChebyshevInterval()
    D=Derivative(d)
    A=D^2-I
    @time κ=nullspace(A)
    @test length(κ) == 2
    c=[κ(0.);κ'(0.)]\[exp(0.),exp(0.)]
    u=(κ*c)[1]
    @test u(1.0) ≈ ℯ
    d=(-50..5.)
    x=Fun(identity,d)
    D=Derivative(d)
    @time u=nullspace(D^2-x)
    c=[u(leftendpoint(d)); u(rightendpoint(d))]\[airyai.(endpoints(d))...]
    @test norm((u*c)[1]-Fun(airyai,d))<10000eps()
end
