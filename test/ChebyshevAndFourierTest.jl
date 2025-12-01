using Test
using ApproxFun

@testset "Chebyshev and Fourier" begin
    @test norm(Fun(x->Fun(cos,Fourier(-π .. π),20)(x),20)-Fun(cos,20)) <100eps()
    @test norm(Fun(x->Fun(cos,Fourier(-π .. π))(x))-Fun(cos)) <100eps()
    @test norm(Fun(x->Fun(cos,Laurent)(x))-Fun(cos)) <100eps()
end
