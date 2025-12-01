using Test
using ApproxFun

@testset "Sampling" begin
    ff=(x,y)->(x-y)^2*exp(-x^2/2-y^2/2)
    f=Fun(ff,Domain(-4..4)^2)
    r=ApproxFun.sample(f,5000)
    g=sum(f,1)/sum(f)
    @test  g(0.1) ≈ 0.2004758624973169
    f = Fun(x -> exp(-x^2/2),-5..5)
    g = cumsum(f)
    @test g(ApproxFun.bisectioninv(g,0.5)) ≈ 0.5
end
