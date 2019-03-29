@testset "Full Singularities Tests" begin
    @time @test norm(Fun(exp,Legendre(0..1))+sqrt(Fun(0..1))) ≈ 2.491141949903508
end

