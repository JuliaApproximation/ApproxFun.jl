using Test
using ApproxFun

@testset "nonlinear pow" begin
    x  = Fun(identity, 0..10)
    n = 2
    u₀ = 0.0x # initial guess
    N = u -> [u'(0), u(0)-1, x*u'' + 2*u' + x*u^n]
    u = newton(N, u₀) # perform Newton iteration in function space
    @test u(0.1) ≈ 0.9983349985461484
    n = 4.5
    u₀ = 0.0x+1 # initial guess
    N = u -> [u'(0), u(0)-1, x*u'' + 2*u' + x*u^n]
    u = newton(N, u₀) # perform Newton iteration in function space
    @test u(0.1) ≈ 0.9983370741307388
end
