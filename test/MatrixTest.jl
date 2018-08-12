using ApproxFun, Test
    import ApproxFun: Block


@testset "RaggedMatrix" begin
    cols=Int[rand(1:k+2) for k=1:5]
    B=ApproxFun.rrand(Float64,maximum(cols),cols)
    cols=Int[rand(1:k+2) for k=1:size(B,1)]
    A=ApproxFun.rrand(Float64,maximum(cols),cols)
    @test full(A)*full(B) ≈ full(A*B)

    @test ApproxFun.RaggedMatrix(B) === B
    @test ApproxFun.RaggedMatrix{Float64}(B) === B
    @test full(ApproxFun.RaggedMatrix{ComplexF64}(B)) == Matrix{ComplexF64}(full(B))

    B = ApproxFun.brand(10,10,2,3)
    @test Matrix(B) == Matrix(ApproxFun.RaggedMatrix(B))
    @test ApproxFun.RaggedMatrix(B) == ApproxFun.RaggedMatrix{Float64}(B)
    @test ApproxFun.RaggedMatrix(ApproxFun.BandedMatrix{ComplexF64}(B)) == ApproxFun.RaggedMatrix{ComplexF64}(B)
end
