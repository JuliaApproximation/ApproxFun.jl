"""
    Represent a symmetric block arrowhead matrix:
    [A   B
     B^⊤ C]
"""
struct SymBlockArrowHead{T,M<:AbstractMatrix{T}} <: AbstractMatrix{T}
    A::M
    B::Matrix{T}
    C::Matrix{T}
    function SymBlockArrowHead(A::M, B::Matrix{T}, C::Matrix{T}) where {T, M<:AbstractMatrix{T}}
        @assert size(A, 1) == size(A, 2) == size(B, 1)
        @assert size(B, 2) == size(C, 1) == size(C, 2)
        new{T, M}(A, B, C)
    end
end

size(A::SymBlockArrowHead) = (s = sum(size(A.B)); (s, s))

function getindex(A::SymBlockArrowHead{T}, i::Integer, j::Integer) where T
    @boundscheck checkbounds(A, i, j)
    if i ≤ size(A.A, 1) && j ≤ size(A.A, 2)
        A.A[i,j]
    elseif i ≤ size(A.A, 1) && j > size(A.A, 2)
        return A.B[i, j - size(A.A, 2)]
    elseif j ≤ size(A.A, 2) && i > size(A.A, 1)
        return A.B[j, i - size(A.A, 1)]
    elseif size(A.A, 1) < i ≤ size(A, 1) && size(A.A, 2) < j ≤ size(A, 2)
        return A.C[i - size(A.A, 1), j - size(A.A, 2)]
    else
        return zero(T)
    end
end

function setindex!(A::SymBlockArrowHead{T}, v, i::Integer, j::Integer) where T
    if i ≤ size(A.A, 1) && j ≤ size(A.A, 2)
        A.A[i,j] = v
    elseif i ≤ size(A.A, 1) && j > size(A.A, 2)
        return A.B[i, j - size(A.A, 2)] = v
    elseif j ≤ size(A.A, 2) && i > size(A.A, 1)
        return A.B[j, i - size(A.A, 1)] = v
    elseif size(A.A, 1) < i ≤ size(A, 1) && size(A.A, 2) < j ≤ size(A, 2)
        return A.C[i - size(A.A, 1), j - size(A.A, 2)] = v
    else
        throw(error("Index out of reach."))
    end
end

function computeHouseholder!(A::SymBlockArrowHead{T,SymBandedPlusBulge{T}}, w::Vector{T}, row::Int, col::Int) where T
    b = A.A.b
    corr = zero(T)
    fill!(w, corr)
    @inbounds for i = max(row-b, 1):row
        corr += abs2(A[i, col])
        w[i] = A[i, col]
    end
    normw = sqrt(corr)
    corr = copysign(normw, A[row, col])
    w[row] += corr
    normw = sqrt(normw^2+corr*muladd(T(2),w[row],-corr))
    if normw == zero(normw)
        return w
    else
        return LinearAlgebra.__normalize!(w, normw)
    end
end

function applyHouseholder!(w::Vector{T}, wQ::Vector{T}, Q::Matrix{T}, bulge::Int, b::Int, row::Int) where T
    sw = max(bulge - b, 1)
    fw = min(bulge, size(Q, 1))
    s = 1
    f = size(Q, 1)

    # wQ = (w^⊤ Q)^⊤
    fill!(wQ, zero(T))
    @inbounds for j = s:f
        wQj = zero(T)
        @simd for i = sw:fw
            wQj += w[i]*Q[i,j]
        end
        wQ[j] = wQj
    end

    # H Q = Q - 2/dot(w, w) (w (w^⊤ Q))
    # twodivnrmw = 2/dot(w, w)
    twodivnrmw = T(2)
    @inbounds for j = s:f
        wQj = wQ[j]
        @simd for i = sw:fw
            Q[i,j] -= twodivnrmw*w[i]*wQj
        end
    end
    #Q .= Q - 2*w*(w'Q)

    Q
end

function similarity!(w::Vector{T}, v::Vector{T}, Aw::Vector{T}, A::SymBlockArrowHead{T,SymBandedPlusBulge{T}}, row::Int) where T
    bulge = A.A.bulge
    b = A.A.b
    sw = max(row-b, 1)
    fw = min(row, size(A, 2))
    sv = max(bulge - 3b, 1)
    fv = min(bulge, size(A, 1))
    sv1 = size(A, 2) - b + 1
    fv1 = size(A, 2)
    s = min(sv, sw)
    f = max(fv, fw)

    # Aw = A*w
    fill!(Aw, zero(T))
    @inbounds for j = sw:fw
        wj = w[j]
        @simd for i = sv:fv
            Aw[i] += A[i,j]*wj
        end
        @simd for i = sv1:fv1
            Aw[i] += A[i,j]*wj
        end
    end

    # v = Aw - dot(w, Aw)/dot(w, w)*w
    fill!(v, zero(T))
    # cst = -dot(w, Aw)/dot(w, w)
    cst = -dot(w, Aw)
    @inbounds @simd for i = sv:fv
        v[i] = muladd(cst, w[i], Aw[i])
    end
    @inbounds @simd for i = sv1:fv1
        v[i] = muladd(cst, w[i], Aw[i])
    end

    # H A H^⊤ = A - 2/dot(w, w) ( vw^⊤ + wv^⊤ )
    # twodivnrmw = 2/dot(w, w)
    twodivnrmw = T(2)
    @inbounds for j = s:f
        vj = v[j]
        wj = w[j]
        @simd for i = s:f
            A[i,j] -= twodivnrmw*(v[i]*wj+w[i]*vj)
        end
    end
    @inbounds for j = sv1:fv1
        vj = v[j]
        wj = w[j]
        @simd for i = s:f
            A[i,j] -= twodivnrmw*(v[i]*wj+w[i]*vj)
        end
    end

    A
end

function symblockarrowhead2symbanded!(w::Vector{T}, v::Vector{T}, Aw::Vector{T}, A::SymBlockArrowHead{T,SymBandedPlusBulge{T}}) where T
    n = size(A, 2)
    b = A.A.b

    for j = 1:b-1
        computeHouseholder!(A, w, b+1-j, n+1-j)
        similarity!(w, v, Aw, A, b+1-j)
    end

    for m = 2:floor(Int, n÷b)-1
        A.A.bulge = m*A.A.b
        for j = 1:b
            computeHouseholder!(A, w, m*b+1-j, n+1-j)
            similarity!(w, v, Aw, A, m*b+1-j)
        end
        chasebulge!(w, v, Aw, A.A)
    end

    A
end

function symblockarrowhead2symbanded!(A::SymBlockArrowHead{T,SymBandedPlusBulge{T}}) where T
    w = zeros(T, size(A, 2))
    v = zero(w)
    Aw = zero(w)
    symblockarrowhead2symbanded!(w, v, Aw, A)
end

function symblockarrowhead2symbanded!(w::Vector{T}, v::Vector{T}, Aw::Vector{T}, A::SymBlockArrowHead{T,SymBandedPlusBulge{T}}, wQ::Vector{T}, Q::Matrix{T}) where T
    n = size(A, 2)
    b = A.A.b

    for j = 1:b-1
        computeHouseholder!(A, w, b+1-j, n+1-j)
        applyHouseholder!(w, wQ, Q, b+1-j, b, n+1-j)
        similarity!(w, v, Aw, A, b+1-j)
    end

    for m = 2:floor(Int, n÷b)-1
        A.A.bulge = m*A.A.b
        for j = 1:b
            computeHouseholder!(A, w, m*b+1-j, n+1-j)
            applyHouseholder!(w, wQ, Q, m*b+1-j, b, n+1-j)
            similarity!(w, v, Aw, A, m*b+1-j)
        end
        chasebulge!(w, v, Aw, A.A, wQ, Q)
    end

    A
end

function symblockarrowhead2symbanded!(A::SymBlockArrowHead{T,SymBandedPlusBulge{T}}, Q::Matrix{T}) where T
    w = zeros(T, size(A, 2))
    v = zero(w)
    Aw = zero(w)
    wQ = zero(w)
    symblockarrowhead2symbanded!(w, v, Aw, A, wQ, Q)
end

function structured_cholesky!(A::SymBlockArrowHead{T,Diagonal{T,V}}) where {T,V}
    D = A.A
    d = D.diag
    d .= sqrt.(d)
    B = A.B
    ldiv!(D, B)
    A.C .= cholesky(A.C-B'B).factors
    UpperTriangular(A)
end

function structured_invcholesky!(A::SymBlockArrowHead{T, Diagonal{T,V}}) where {T,V}
    D = A.A
    d = D.diag
    d .= sqrt.(d)
    B = A.B
    ldiv!(D, B)
    A.C .= cholesky(A.C-B'B).factors
    # A, B, C are the old D^{1/2}, E, F.
    d .= inv.(d)
    A.C .= inv(UpperTriangular(A.C)).data
    B .= -lmul!(D, B*A.C)
    UpperTriangular(A)
end

function congruence!(L::SymBlockArrowHead{T,M}, C::SymBlockArrowHead{T, Diagonal{T,V}}) where {T, M, V}
    R = structured_invcholesky!(C)
    A = L.A
    B = L.B
    CC = L.C
    D = R.data.A
    E = R.data.B
    F = UpperTriangular(R.data.C)
    DAD = D*A*D
    DAEpDBF = D*(A*E + B*F)
    EtBF = E'B*F
    EtAEpEtBFpFtBtEpFtCF = E'A*E + F'CC*F + (EtBF+EtBF')
    A .= DAD
    B .= DAEpDBF
    CC .= EtAEpEtBFpFtBtEpFtCF
    L, C
end
