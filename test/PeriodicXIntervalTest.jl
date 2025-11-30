using Test
using ApproxFun
using Random
using LinearAlgebra
using ApproxFunBase

@testset "periodic x interval" begin
    dθ=PeriodicSegment(-2.,2.)
    dt=Interval(0,1.)
    d=dθ×dt
    Dθ=Derivative(d,[1,0])
    Dt=Derivative(d,[0,1])
    u0=Fun(θ->exp(-20θ^2), dθ)
    ε = 0.1
    u=\([I⊗ldirichlet(dt); Dt-ε*Dθ^2-Dθ], [u0; 0.]; tolerance=1E-4)
    @test ≈(u(0.1,0.2),0.3103472600253807;atol=1E-2)
    A=Dt+Dθ
    ApproxFunBase.TestUtils.testbandedblockbandedoperator(A)
    u=\([I⊗ldirichlet(dt); Dt+Dθ], [u0; 0.0]; tolerance=1E-6)
    @test ≈(u(0.2,0.1),u0(0.1);atol=1E-6)

    d=PeriodicSegment() × ChebyshevInterval()
    u_ex=Fun((x,y)->real(cos(x+im*y)),d)
    @test u_ex(1.0,0.1) ≈ real(cos(1.0+im*0.1)) atol=10eps()
    B=Dirichlet(Space(d))
    @test B.order == 0
    g=Fun((x,y)->real(cos(x+im*y)),rangespace(B))
    @test norm((B*u_ex-g).coefficients) < 100eps()
    testbandedblockbandedoperator(Laplacian(d))
    @time u=[B;Laplacian(d)]\[g;0.]
    @test u(.1,.2) ≈ real(cos(.1+.2im))
end
