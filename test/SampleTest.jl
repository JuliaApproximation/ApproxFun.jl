using Test
using ApproxFun

@testset "sample" begin
    f=Fun(exp)
    sample(f,100000)
    @time sample(f,100000)
end
