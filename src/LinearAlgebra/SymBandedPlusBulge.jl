"""
    Represent a symmetric banded matrix plus a bulge:
    [ □ ◺
      ◹ □ ◺
        ◹ □ □ ◺
          □ □ □
          ◹ □ □ ◺
              ◹ □ ◺
                ◹ □ ]
"""
mutable struct SymBandedPlusBulge{T,M<:BandedMatrix{T}} <: AbstractMatrix{T}
    A::M
    b::Int
    bulge::Int
end

SymBandedPlusBulge(A::AbstractMatrix, b::Int, bulge::Int) = SymBandedPlusBulge(BandedMatrix(A, (2b, 2b)), b, bulge)

size(A::SymBandedPlusBulge) = size(A.A)

function getindex(A::SymBandedPlusBulge{T}, i::Integer, j::Integer) where T
    b = A.b
    bulge = A.bulge
    col1 = bulge - 2b
    col2 = bulge - b
    row1 = col1 - b
    row2 = col2 - 2b
    AA = A.A
    if -b ≤ i-j ≤ b
        if j ≥ i
            return AA[i,j]
        else
            return AA[j,i]
        end
    elseif 1 ≤ j-col1 ≤ b && 1 ≤ i-row1 ≤ b && j-col1 > i-row1
        return AA[i,j]
    elseif 1 ≤ j-col2 ≤ b && 1 ≤ i-row2-(j-col2-1) ≤ b
        return AA[i,j]
    elseif 1 ≤ i-col1 ≤ b && 1 ≤ j-row1 ≤ b && i-col1 > j-row1
        return AA[j,i]
    elseif 1 ≤ i-col2 ≤ b && 1 ≤ j-row2-(i-col2-1) ≤ b
        return AA[j,i]
    else
        return zero(T)
    end
end

setindex!(A::SymBandedPlusBulge{T}, v, i::Integer, j::Integer) where T = (b = A.b; -2b ≤ i-j ≤ 2b && setindex!(A.A, v, i, j))

function computeHouseholder!(A::SymBandedPlusBulge{T}, w::Vector{T}, col::Int) where T
    b = A.b
    corr = zero(T)
    fill!(w, corr)
    row = col - b
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

function applyHouseholder!(w::Vector{T}, wQ::Vector{T}, Q::Matrix{T}, bulge::Int, b::Int) where T
    sw = max(bulge - 2b, 1)
    fw = min(bulge - b, size(Q, 1))
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

function similarity!(w::Vector{T}, v::Vector{T}, Aw::Vector{T}, A::SymBandedPlusBulge{T}) where T
    bulge = A.bulge
    b = A.b
    sw = max(bulge - 2b, 1)
    fw = min(bulge - b, size(A, 1))
    sv = max(bulge - 3b, 1)
    fv = min(bulge, size(A, 1))
    s = min(sv, sw)
    f = max(fv, fw)

    # Aw = A*w
    fill!(Aw, zero(T))
    @inbounds for j = sw:fw
        wj = w[j]
        @simd for i = sv:fv
            Aw[i] += A[i,j]*wj
        end
    end

    # v = Aw - dot(w, Aw)/dot(w, w)*w
    fill!(v, zero(T))
    # cst = -dot(w, Aw)/dot(w, w)
    cst = -dot(w, Aw)
    @inbounds for i = sv:fv
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

    A.bulge -= 1
    A
end

function chasebulge!(w::Vector{T}, v::Vector{T}, Aw::Vector{T}, A::SymBandedPlusBulge{T}) where T
    bulge = A.bulge
    b = A.b
    for k = bulge:-1:b+1
        computeHouseholder!(A, w, k)
        similarity!(w, v, Aw, A)
    end

    A
end

function chasebulge!(A::SymBandedPlusBulge{T}) where T
    w = zeros(T, size(A, 2))
    v = zero(w)
    Aw = zero(w)
    chasebulge!(w, v, Aw, A)
end

function chasebulge!(w::Vector{T}, v::Vector{T}, Aw::Vector{T}, A::SymBandedPlusBulge{T}, wQ::Vector{T}, Q::Matrix{T}) where T
    bulge = A.bulge
    b = A.b
    for k = bulge:-1:b+1
        computeHouseholder!(A, w, k)
        applyHouseholder!(w, wQ, Q, k, b)
        similarity!(w, v, Aw, A)
    end

    A, Q
end

function chasebulge!(A::SymBandedPlusBulge{T}, Q::Matrix{T}) where T
    w = zeros(T, size(A, 2))
    v = zero(w)
    Aw = zero(w)
    wQ = zero(w)
    chasebulge!(w, v, Aw, A, wQ, Q)
end
