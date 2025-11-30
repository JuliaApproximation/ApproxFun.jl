using Test
using ApproxFun

@testset "Laplace in a strip" begin
    d = PeriodicSegment() × ChebyshevInterval()
    g=Fun((x,y)->real(cos(x+im*y)),∂(d))
    @test g(0.1,1.0) ≈ real(cos(0.1+im))
    @test g(0.1,-1.0) ≈ real(cos(0.1-im))
    v=[g;0]
    @test v(0.1,-1) ≈ [real(cos(0.1-im));0]
    A=[Dirichlet(d);Laplacian(d)]
    a = space(v)
    b = rangespace(A)
    @test Fun(component(v[1],1), component(b[1],1))(0.1,-1.0) ≈ v(0.1,-1.0)[1]
    @test Fun(component(v[1],2), component(b[1],2))(0.1,1.0) ≈ v(0.1,1.0)[1]
    @test ApproxFun.default_Fun(v[1] , b[1])(0.1,1.0) ≈ v(0.1,1.0)[1]
end
