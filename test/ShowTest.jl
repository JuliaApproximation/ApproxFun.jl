using Test
using ApproxFun

@testset "show" begin
    op = Derivative(Chebyshev())
    io = IOBuffer()
    @test summary(io, op) isa Nothing
    @test contains(String(take!(io)), " : $(domainspace(op)) → $(rangespace(op))")
    show(io, op)
    @test contains(String(take!(io)), " : $(domainspace(op)) → $(rangespace(op))")
    @test summary(io, ApproxFun.ArraySpace(Chebyshev(), 2)) isa Nothing
    @test contains(String(take!(io)), "ArraySpace")
    Q = QuotientSpace(Dirichlet(Chebyshev()))
    @test startswith(repr(Q), "Chebyshev() /")
    show(io, MIME"text/plain"(), Q)
    s = String(take!(io))
    @test startswith(s, "Chebyshev() /")
    f = Fun(Chebyshev()^2, [1,3,4])
    @test contains(repr(f), repr(space(f)))
    @test contains(repr(f), repr(coefficients(f)))
end
