using Test
using ApproxFun

@testset "Mix Fourier-Chebyshev (#602)" begin
    s = Chebyshev(-π..π)
    a = Fun(t-> 1+sin(cos(2t)), s)
    L = Derivative() + a
    f = Fun(t->exp(sin(10t)), s)
    B = periodic(s,0)
    @time uChebyshev = [B;L] \ [0.;f]
    s = Fourier(-π..π)
    a = Fun(t-> 1+sin(cos(2t)), s)
    L = Derivative() + a
    f = Fun(t->exp(sin(10t)), s)
    @time uFourier = L\f
    @test norm(uFourier-uChebyshev) ≤ 100eps()
end
