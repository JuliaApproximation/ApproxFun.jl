using ApproxFun, BandedMatrices, SpecialFunctions, Test
    import ApproxFun: bandwidth

#
# Fix for BigFloat stackoverflow
# https://github.com/JuliaLang/julia/pull/30952
#
Base._range(a::T, step::T, ::Nothing, len::Integer) where {T <: AbstractFloat} =
	    Base._rangestyle(Base.OrderStyle(T), Base.ArithmeticStyle(T), a, step, len)

@testset "Eigenvalue problems" begin
    @testset "Negative Laplacian with Neumann boundary conditions" begin
        #
        # -𝒟² u = λu,  u'(±1) = 0.
        #
        d = Segment(-1..1)
        S = Ultraspherical(0.5, d)
        NS = NormalizedPolynomialSpace(S)
        L = -Derivative(S, 2)
        C = Conversion(domainspace(L), rangespace(L))
        B = Neumann(S)
        QS = QuotientSpace(B)
        Q = Conversion(QS, S)
        D1 = Conversion(S, NS)
        D2 = Conversion(NS, S)
        R = D1*Q
        P = cache(PartialInverseOperator(C, (0, bandwidth(L, 1) + bandwidth(R, 1) + bandwidth(C, 2))))
        A = R'D1*P*L*D2*R
        B = R'R

        # Currently a hack to avoid a LAPACK calling bug with a pencil (A, B)
        # with smaller bandwidth in A.
        n = 50
        SA = Symmetric(Matrix(A[1:n,1:n]), :L)
        SB = Symmetric(Matrix(B[1:n,1:n]), :L)
        λ = eigvals(SA, SB)

        @test λ[1:round(Int, 2n/5)] ≈ (π^2/4).*(0:round(Int, 2n/5)-1).^2
    end

    @testset "Schrödinger with piecewise-linear potential with Dirichlet boundary conditions" begin
        #
        # [-𝒟² + V] u = λu,  u(±1) = 0,
        #
        # where V = 100|x|.
        #
        d = Segment(-1..0)∪Segment(0..1)
        S = PiecewiseSpace(Ultraspherical.(0.5, d.domains))
        NS = PiecewiseSpace(NormalizedUltraspherical.(0.5, d.domains))
        V = 100Fun(abs, S)
        L = -Derivative(S, 2) + V
        C = Conversion(domainspace(L), rangespace(L))
        B = [Dirichlet(S); continuity(S, 0:1)]
        QS = QuotientSpace(B)
        Q = Conversion(QS, S)
        D1 = Conversion(S, NS)
        D2 = Conversion(NS, S)
        R = D1*Q
        P = cache(PartialInverseOperator(C, (0, bandwidth(L, 1) + bandwidth(R, 1) + bandwidth(C, 2))))
        A = R'D1*P*L*D2*R
        B = R'R

        n = 100
        SA = Symmetric(A[1:n,1:n], :L)
        SB = Symmetric(B[1:n,1:n], :L)
        λ = eigvals(SA, SB)

        @test λ[1] ≈ parse(BigFloat, "2.19503852085715200848808942880214615154684642693583513254593767079468401198338e+01")
    end

    @testset "Schrödinger with piecewise-constant potential with Dirichlet boundary conditions" begin
        #
        # [-𝒟² + V] u = λu,  u(±1) = 0,
        #
        # where V = 1000[χ_[-1,-1/2](x) + χ_[1/2,1](x)].
        #
        d = Segment(-1..(-0.5))∪Segment(-0.5..0.5)∪Segment(0.5..1)
        S = PiecewiseSpace(Ultraspherical.(0.5, d.domains))
        NS = PiecewiseSpace(NormalizedUltraspherical.(0.5, d.domains))
        V = Fun(x->abs(x) ≥ 1/2 ? 1000 : 0, S)
        L = -Derivative(S, 2) + V
        C = Conversion(domainspace(L), rangespace(L))
        B = [Dirichlet(S); continuity(S, 0:1)]
        QS = QuotientSpace(B)
        Q = Conversion(QS, S)
        D1 = Conversion(S, NS)
        D2 = Conversion(NS, S)
        R = D1*Q
        P = cache(PartialInverseOperator(C, (0, bandwidth(L, 1) + bandwidth(R, 1) + bandwidth(C, 2))))
        A = R'D1*P*L*D2*R
        B = R'R

        n = 150
        SA = Symmetric(A[1:n,1:n], :L)
        SB = Symmetric(B[1:n,1:n], :L)
        λ = eigvals(SA, SB)
        # From Lee--Greengard (1997).
        λtrue = [2.95446;5.90736;8.85702;11.80147].^2
        @test norm((λ[1:4] - λtrue)./λ[1:4]) < 1e-5
    end

    @testset "BigFloat negative Laplacian with Dirichlet boundary conditions" begin
        #
        # -𝒟² u = λu,  u(±1) = 0.
        #
        d = Segment(big(-1.0)..big(1.0))
        S = Ultraspherical(big(0.5), d)
        NS = NormalizedPolynomialSpace(S)
        L = -Derivative(S, 2)
        C = Conversion(domainspace(L), rangespace(L))
        B = Dirichlet(S)
        QS = QuotientSpace(B)
        Q = Conversion(QS, S)
        D1 = Conversion(S, NS)
        D2 = Conversion(NS, S)
        R = D1*Q
        P = cache(PartialInverseOperator(C, (0, bandwidth(L, 1) + bandwidth(R, 1) + bandwidth(C, 2))))
        A = R'D1*P*L*D2*R
        B = R'R

        n = 300
        SA = Symmetric(A[1:n,1:n], :L)
        SB = Symmetric(B[1:n,1:n], :L)
        BSA = BandedMatrix(SA)
        BSB = BandedMatrix(SB)
        begin
            v = zeros(BigFloat, n)
            v[1] = 1
            v[3] = 1/256
            λ = big(0.0)
            for _ = 1:7
                λ = dot(v, BSA*v)/dot(v, BSB*v)
                OP = Symmetric(SA.data-λ*SB.data, :L)
                F = ldlt!(OP)
                ldiv!(F, v)
                normalize!(v)
            end
            @test λ ≈ big(π)^2/4
        end
    end

    @testset "Skew differentiation matrices" begin
        #
        # 𝒟 u = λ u,  u(1) + u(-1) = 0.
        #
        d = Segment(-1..1)
        S = Ultraspherical(0.5, d)
        NS = NormalizedPolynomialSpace(S)
        Lsk = Derivative(S)
        Csk = Conversion(domainspace(Lsk), rangespace(Lsk))
        B = Evaluation(S, -1) + Evaluation(S, 1)
        QS = PathologicalQuotientSpace(B)
        Q = Conversion(QS, S)
        D1 = Conversion(S, NS)
        D2 = Conversion(NS, S)
        R = D1*Q
        Psk = cache(PartialInverseOperator(Csk, (0, 4 + bandwidth(Lsk, 1) + bandwidth(R, 1) + bandwidth(Csk, 2))))
        Ask = R'D1*Psk*Lsk*D2*R
        Bsk = R'R

        n = 100
        λim = imag(sort!(eigvals(Matrix(tril(Ask[1:n,1:n], 1)), Matrix(tril(Bsk[1:n,1:n], 2))), by = abs))

        @test abs.(λim[1:2:round(Int, 2n/5)]) ≈ π.*(0.5:round(Int, 2n/5)/2)
        @test abs.(λim[2:2:round(Int, 2n/5)]) ≈ π.*(0.5:round(Int, 2n/5)/2)
    end

    @testset "Complex spectra of principal finite sections of an ultraspherical discretization of a self-adjoint problem" begin
        #
        # -𝒟² u = λ u,  u'(±1) = 0.
        #
        S = Chebyshev()
        L = -Derivative(S, 2)
        B = Neumann(S)
        C = Conversion(domainspace(L), rangespace(L))
        n = 10 # ≥ 10 appears to do the trick
        λ = eigvals(Matrix([B;L][1:n,1:n]), Matrix([B-B;C][1:n,1:n]))
        @test eltype(λ) == Complex{Float64}
    end
end
