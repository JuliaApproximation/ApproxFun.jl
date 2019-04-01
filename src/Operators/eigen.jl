const sqrt2m1 = sqrt(2)-1

function symmetric_eigen(L::Operator{T}, n::Int; verbose::Bool = false, tol::Real = eps(T), CONST::Real = sqrt2m1) where T
    # ndef is the number of deflated eigenvalues
    # ALLEIGS is a container for the converged eigenvalues
    # b is the bandwidth of the symmetric banded operator
    # (k, k) are the dimensions of the deflation window
    ds = domainspace(L)
    rs = rangespace(L)
    ndef = 0
    ALLEIGS = T[]
    ALLQ = EmbeddedBlockOperator{T,Matrix{T},typeof(ds),typeof(rs)}[]
    b = maximum(bandwidths(L))
    k = 2b

    # We start by determining a good deflation window.

    Lk = Symmetric(L[1:k,1:k])
    Lb = L[k+1:k+b,1:k]
    Λ, Q = eigen(Lk)
    M = Lb*Q

    j = 0
    while j < 1
        k *= 2
        Lk = Symmetric(L[1:k,1:k])
        Lb = L[k+1:k+b,1:k]
        Λ, Q = eigen(Lk)
        M = Lb*Q
        j = deflationcheck(M, Λ, tol^2, verbose)÷b*b
        verbose && println("This is the size of the deflation window: ",k," and this is the number of deflated eigenvalues: ",j)
    end

    SBA = SymBlockArrowHead(SymBandedPlusBulge(Diagonal(Λ[j+1:k]), b, b), Matrix(M[:,j+1:k]'), Matrix(L[k+1:k+b,k+1:k+b]))

    push!(ALLEIGS, Λ[1:j]...)
    push!(ALLQ, EmbeddedBlockOperator(Q, ndef, ds, rs))
    ndef += j

    w = zeros(T, k+b)
    v = zero(w)
    Lw = zero(w)
    wQ = zero(w)

    while ndef < n
        # With the dimensions of the deflation window, (k, k), found & fixed,
        # we now create a SymBlockArrowHead matrix upon which we are able to
        # perform similarity transformations to return to banded form.

        Q = Matrix{T}(I, k, k)
        symblockarrowhead2symbanded!(w, v, Lw, SBA, wQ, Q)
        push!(ALLQ, EmbeddedBlockOperator(Matrix(Q'), ndef, ds, rs))

        replenishLk!(Lk, SBA, L, ndef)
        replenishLb!(Lb, L, ndef, k, b)

        if j < CONST*k
            verbose && println("Failed to deflate! Doubling the dimensions of the deflation window.")
            Lk = [Lk L[ndef+1:ndef+k,ndef+k+1:ndef+2k]; L[ndef+k+1:ndef+2k,ndef+1:ndef+2k]]
            k *= 2
            Lk = SymBandedMatrix(Lk, b, 1:k)
            Lb = L[ndef+k+1:ndef+k+b, ndef+1:ndef+k]
            M = zeros(T, b, k)
            w = zeros(T, k+b)
            v = zero(w)
            Lw = zero(w)
            wQ = zero(w)
        end

        Λ, Q = eigen(Lk)

        mul!(fill!(M, zero(T)), Lb, Q)
        j = deflationcheck(M, Λ, tol^2, verbose)÷b*b
        verbose && println("This is the size of the deflation window: ",k," and this is the number of deflated eigenvalues: ",j)

        SBA = SymBlockArrowHead(SymBandedPlusBulge(Diagonal(Λ[j+1:k]), b, b), Matrix(M[:,j+1:k]'), Matrix(L[ndef+k+1:ndef+k+b,ndef+k+1:ndef+k+b]))

        push!(ALLEIGS, Λ[1:j]...)
        push!(ALLQ, EmbeddedBlockOperator(Q, ndef, ds, rs))
        ndef += j
    end

    ALLEIGS, ALLQ
end

function symmetric_eigen(L::Operator{T}, C::Operator{T}, n::Int; verbose::Bool = false, tol::Real = eps(T), CONST::Real = sqrt2m1) where T
    # ndef is the number of deflated eigenvalues
    # ALLEIGS is a container for the converged eigenvalues
    # b is the bandwidth of the symmetric banded operator
    # (k, k) are the dimensions of the deflation window
    ds = domainspace(L)
    rs = rangespace(L)
    ndef = 0
    ALLEIGS = T[]
    ALLV = EmbeddedBlockOperator{T,Matrix{T},typeof(ds),typeof(rs)}[]
    b = max(maximum(bandwidths(L)), maximum(bandwidths(C)))
    k = 2b

    # We start by determining a good deflation window.

    Lk = Symmetric(L[1:k,1:k])
    Lb = L[k+1:k+b,1:k]
    Ck = Symmetric(C[1:k,1:k])
    Cb = C[k+1:k+b,1:k]
    Λ, V = eigen(Lk, Ck)
    M = Lb*V
    N = Cb*V

    j = 0
    while j < 2b
        k *= 2
        Lk = Symmetric(L[1:k,1:k])
        Lb = L[k+1:k+b,1:k]
        Ck = Symmetric(C[1:k,1:k])
        Cb = C[k+1:k+b,1:k]
        Λ, V = eigen(Lk, Ck)
        M = Lb*V
        N = Cb*V
        j = deflationcheck(M, N, Λ, tol^2, verbose)÷b*b
        verbose && println("This is the size of the deflation window: ",k," and this is the number of deflated eigenvalues: ",j)
    end

    SBAL = SymBlockArrowHead(SymBandedPlusBulge(Diagonal(Λ[j+1:k]), b, b), Matrix(M[:,j+1:k]'), Matrix(L[k+1:k+b,k+1:k+b]))
    SBAC = SymBlockArrowHead(Diagonal(ones(T, k-j)), Matrix(N[:,j+1:k]'), Matrix(C[k+1:k+b,k+1:k+b]))

    push!(ALLEIGS, Λ[1:j]...)
    push!(ALLV, EmbeddedBlockOperator(V, ndef, ds, rs))
    ndef += j

    w = zeros(T, k+b)
    v = zero(w)
    Lw = zero(w)
    wQ = zero(w)

    while ndef < n
        # With the dimensions of the deflation window, (k, k), found & fixed,
        # we now create a SymBlockArrowHead matrix upon which we are able to
        # perform similarity transformations to return to banded form.

        congruence!(SBAL, SBAC)
        push!(ALLV, EmbeddedBlockOperator(Matrix(UpperTriangular(SBAC)), ndef, ds, rs))

        Q = Matrix{T}(I, k, k)
        symblockarrowhead2symbanded!(w, v, Lw, SBAL, wQ, Q)
        push!(ALLV, EmbeddedBlockOperator(Matrix(Q'), ndef, ds, rs))

        replenishLk!(Lk, SBAL, SBAC, L, ndef)
        replenishLb!(Lb, L, ndef, k, b)
        replenishCk!(Ck, SBAC, C, b, ndef)
        replenishLb!(Cb, C, ndef, k, b)

        if j < CONST*k
            verbose && println("Failed to deflate! Doubling the dimensions of the deflation window.")

            Lk = [Lk L[ndef+1:ndef+k,ndef+k+1:ndef+2k]; L[ndef+k+1:ndef+2k,ndef+1:ndef+2k]]
            Ck = [Ck C[ndef+1:ndef+k,ndef+k+1:ndef+2k]; C[ndef+k+1:ndef+2k,ndef+1:ndef+2k]]
            k *= 2
            Lk = SymBandedMatrix(Lk, b, 1:k)
            Ck = SymBandedMatrix(Ck, b, 1:k)
            Lb = L[ndef+k+1:ndef+k+b, ndef+1:ndef+k]
            Cb = C[ndef+k+1:ndef+k+b, ndef+1:ndef+k]
            M = zeros(T, b, k)
            N = zeros(T, b, k)
            w = zeros(T, k+b)
            v = zero(w)
            Lw = zero(w)
            wQ = zero(w)
        end

        Λ, V = eigen(Lk, Ck)
        mul!(fill!(M, zero(T)), Lb, V)
        mul!(fill!(N, zero(T)), Cb, V)
        j = deflationcheck(M, N, Λ, tol^2, verbose)÷b*b
        verbose && println("This is the size of the deflation window: ",k," and this is the number of deflated eigenvalues: ",j)

        SBAL = SymBlockArrowHead(SymBandedPlusBulge(Diagonal(Λ[j+1:k]), b, b), Matrix(M[:,j+1:k]'), Matrix(L[ndef+k+1:ndef+k+b,ndef+k+1:ndef+k+b]))
        SBAC = SymBlockArrowHead(Diagonal(ones(T, k-j)), Matrix(N[:,j+1:k]'), Matrix(C[ndef+k+1:ndef+k+b,ndef+k+1:ndef+k+b]))

        push!(ALLEIGS, Λ[1:j]...)
        push!(ALLV, EmbeddedBlockOperator(V, ndef, ds, rs))
        ndef += j
    end

    ALLEIGS, ALLV
end



symmetric_eigvals(L::Operator{T}, n::Int; kwargs...) where T = symmetric_eigen(L, n; kwargs...)[1][1:n]
symmetric_eigvals(L::Operator{T}, C::Operator{T}, n::Int; kwargs...) where T = symmetric_eigen(L, C, n; kwargs...)[1][1:n]





# Miscellaneous

function SymBandedMatrix(M::AbstractMatrix{T}, b::Int, ir::UnitRange) where T
    B = BandedMatrix{T}(undef, length(ir), length(ir), b, b)
    shft = first(ir)-1
    @inbounds for i in 1:length(ir)
        for j in colrange(B, i)
            B[i,j] = M[i+shft, j+shft]
        end
    end
    Symmetric(B)
end

function replenishLk!(Lk::Symmetric{T,MB}, SBA::SymBlockArrowHead, L::Operator, ndef::Int) where {T,MB<:BandedMatrix{T}}
    k = size(Lk, 1)
    l = size(SBA.A, 1)
    @inbounds for i = 1:l
        for j in colrange(Lk, i)
            Lk.data[i,j] = SBA[i,j]
        end
    end
    @inbounds for i = l+1:k
        for j in colrange(Lk, i)
            Lk.data[i,j] = L[i+ndef,j+ndef]
        end
    end
    Lk
end

function replenishLk!(A::Symmetric{T,MB}, SBAL::SymBlockArrowHead, SBAC::SymBlockArrowHead, M::Operator, ndef::Int) where {T,MB<:BandedMatrix{T}}
    b = SBAL.A.b
    m, n = size(SBAL.A)
    for i = size(SBAL.A, 1)+1:size(A, 1)
        for j in colrange(A, i)
            A.data[i,j] = M[i+ndef,j+ndef]
        end
    end
    @inbounds for i = size(SBAL, 1)-2b+1:size(SBAL, 1)-b
        for j in colrange(A, i)
            A.data[i,j] = SBAL[i,j]
        end
    end
    @inbounds for i = 1:size(SBAL.A, 1)
        for j in colrange(A, i)
            A.data[i,j] = SBAL[i,j]
        end
    end
    @inbounds for i = size(SBAL.A, 1)+1:size(SBAL, 1)
        for j = i:size(SBAL, 2)
            A.data[i,j] = SBAL[i,j]
        end
    end
    k = size(SBAC.A, 1)
    temp1 = UpperTriangular(SBAC.C)
    temp2 = UpperTriangular(Matrix(M[ndef+k+b+1:ndef+k+2b,ndef+k+1:ndef+k+b]))
    temp3 = temp2*temp1
    @inbounds for i = k+b+1:k+2b
        for j = k+1:k+b
            A.data[i,j] = A.data[j,i] = temp3[i-k-b,j-k]
        end
    end
    A
end

function replenishCk!(A::Symmetric{T,MB}, SBA::SymBlockArrowHead, M::Operator, b::Int, ndef::Int)  where {T,MB<:BandedMatrix{T}}
    m, n = size(SBA.A)
    @inbounds for i = size(SBA.A, 1)+1:size(A, 1)
        for j in colrange(A, i)
            A.data[i,j] = M[i+ndef,j+ndef]
        end
    end
    @inbounds for i = 1:size(SBA.A, 1)+b
        for j in colrange(A, i)
            A.data[i,j] = 0
        end
        A.data[i,i] = 1
    end
    k = size(SBA.A, 1)
    temp1 = UpperTriangular(SBA.C)
    temp2 = UpperTriangular(Matrix(M[ndef+k+b+1:ndef+k+2b,ndef+k+1:ndef+k+b]))
    temp3 = temp2*temp1
    @inbounds for i = k+b+1:k+2b
        for j = k+1:k+b
            A.data[i,j] = A.data[j,i] = temp3[i-k-b,j-k]
        end
    end
    A
end

function replenishLb!(Lb::BandedMatrix{T}, L::Operator{T}, ndef::Int, k::Int, b::Int) where T
    ndefpk = ndef+k
    @inbounds for j = k-b+1:k
        jpndef = j+ndef
        for i in 1:(b+j-k)
            Lb[i,j] = L[i+ndefpk, jpndef]
        end
    end
    Lb
end
