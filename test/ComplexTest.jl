using ApproxFun, Compat.Test

@testset "Complex" begin
    ## Diff
    f=Fun(x->exp(im*x))
        ncoefficients(f)

    @test norm(f'-im*f) < ncoefficients(f)*100*eps()


    @test norm(integrate(f)+im*f-f.coefficients[1]*im) < 100eps()


    @test norm(real(f)-Fun(cos)) < 10eps()

    @test norm(real(f-Fun(cos))) < 10eps()


    ##Check other real domains


    f=Fun(x->exp(im.*x),1..2)
    @test norm(f'-im*f) < 1000eps()
    @test norm(integrate(f)+im*f-f.coefficients[1]*im) < 100eps()
    @test norm(real(f)-Fun(cos,domain(f))) < eps()
    @test norm(real(f-Fun(cos,domain(f)))) < eps()


    ##Check complex domains


    f=Fun(x->exp(im.*x),1im..(2+.5im))

    @test sum(f) ≈ (0.5515167681675808 + 0.6202852564797062im)
    @test f(1im) ≈ exp(im.*im)
    @test f(1+.75im) ≈ exp(im.*(1+.75im))
    @test norm(f'-im*f) < 1000eps()
    @test norm(integrate(f)+im*f-f.coefficients[1]*im) < 100eps()
end
