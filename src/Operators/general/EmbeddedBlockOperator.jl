export EmbeddedBlockOperator

"""
    Represent a square block operator embedded in the identity:
    [ I
        X
          I ]
    where the first identity is k x k and the second is infinite-dimensional.
"""
struct EmbeddedBlockOperator{T,M<:AbstractMatrix{T},DS,RS} <: Operator{T}
    X::M
    temp::Vector{T}
    k::Int
    domainspace::DS
    rangespace::RS
    function EmbeddedBlockOperator{T,M,DS,RS}(X::M, k::Int, ds::DS, rs::RS) where {T,M,DS,RS}
        m, n = size(X)
        @assert m == n
        new{T,M,DS,RS}(X, zeros(T, m), k, ds, rs)
    end
end

EmbeddedBlockOperator(X::M, k::Int, ds, rs) where M = EmbeddedBlockOperator{eltype(X), M, typeof(ds), typeof(rs)}(X, k, ds, rs)
EmbeddedBlockOperator(X::M, k::Int) where M = EmbeddedBlockOperator(X, k, SequenceSpace(), SequenceSpace())

bandwidths(A::EmbeddedBlockOperator) = (size(A.X, 1)-1, size(A.X, 2)-1)
domainspace(A::EmbeddedBlockOperator) = A.domainspace
rangespace(A::EmbeddedBlockOperator) = A.rangespace

function getindex(A::EmbeddedBlockOperator{T}, i::Integer, j::Integer) where T
    k = A.k
    n = size(A.X, 1)
    if i ≤ k
        if j ≤ k
            if i == j
                return one(T)
            else
                return zero(T)
            end
        else
            return zero(T)
        end
    elseif k < i ≤ k+n
        if k < j ≤ k+n
            return A.X[i-k,j-k]
        else
            return zero(T)
        end
    else
        if i == j
            return one(T)
        else
            return zero(T)
        end
    end
end

function lmul!(A::EmbeddedBlockOperator{T}, x::AbstractVector{T}) where T
    X = A.X
    y = A.temp
    k = A.k
    n = size(X, 1)
    @assert length(x) ≥ k + n
    fill!(y, zero(T))
    @inbounds for j = 1:n
        xj = x[k+j]
        for i = 1:n
            y[i] += X[i,j]*xj
        end
    end
    @inbounds for i = k+1:k+n x[i] = y[i-k] end

    x
end

function lmul!(A::AdjointOperator{T,M}, x::AbstractVector{T}) where {T,M<:EmbeddedBlockOperator{T}}
    X = A.op.X
    y = A.op.temp
    k = A.op.k
    n = size(X, 1)
    @assert length(x) ≥ k + n
    fill!(y, zero(T))
    @inbounds for i = 1:n
        yi = zero(eltype(y))
        for j = 1:n
            yi += conj(X[j,i])*x[k+j]
        end
        y[i] += yi
    end
    @inbounds for i = k+1:k+n x[i] = y[i-k] end

    x
end

function lmul!(A::TransposeOperator{T,M}, x::AbstractVector{T}) where {T,M<:EmbeddedBlockOperator{T}}
    X = A.op.X
    y = A.op.temp
    k = A.op.k
    n = size(X, 1)
    @assert length(x) ≥ k + n
    fill!(y, zero(T))
    @inbounds for i = 1:n
        yi = zero(eltype(y))
        for j = 1:n
            yi += X[j,i]*x[k+j]
        end
        y[i] += yi
    end
    @inbounds for i = k+1:k+n x[i] = y[i-k] end

    x
end


function lmul!(A::Vector{M}, x::Vector{T}) where {T,M<:EmbeddedBlockOperator{T}}
    for N = length(A):-1:1
        lmul!(A[N], x)
    end
    x
end

function lmul!(A::Adjoint{AdjointOperator{T,M},Vector{M}}, x::Vector{T}) where {T,M<:EmbeddedBlockOperator{T}}
    for N = 1:length(A)
        lmul!(A[N], x)
    end
    x
end

function lmul!(A::Transpose{TransposeOperator{T,M},Vector{M}}, x::Vector{T}) where {T,M<:EmbeddedBlockOperator{T}}
    for N = 1:length(A)
        lmul!(A[N], x)
    end
    x
end
