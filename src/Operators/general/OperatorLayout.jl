

export HermitianOperator, SymmetricOperator, AdjointOperator, TransposeOperator


struct HermitianOperator{T<:Number,B<:Operator} <: Operator{T}
    op::B
    uplo::Char
end

HermitianOperator(B::Operator{T}, uplo::Symbol=:U) where {T<:Number} = HermitianOperator{T,typeof(B)}(B, char_uplo(uplo))
convert(::Type{Operator{T}},A::HermitianOperator) where {T}=HermitianOperator(convert(Operator{T},A.op), A.uplo)

domainspace(P::HermitianOperator)=domainspace(P.op)
rangespace(P::HermitianOperator)=rangespace(P.op)
domain(P::HermitianOperator)=domain(P.op)
bandwidths(P::HermitianOperator) = (b = bandwidth(P.op, P.uplo == 'L' ? 1 : P.uplo == 'U' ? 2 : 0); (b, b))

function getindex(P::HermitianOperator,k::Integer,j::Integer)
    if P.uplo == 'L'
        if j > k
            conj(P.op[j,k])
        else
            P.op[k,j]
        end
    elseif P.uplo == 'U'
        if j < k
            conj(P.op[j,k])
        else
            P.op[k,j]
        end
    end
end

Hermitian(A::Operator, uplo::Symbol)=HermitianOperator(A, uplo)
Hermitian(A::Operator)=HermitianOperator(A)


struct SymmetricOperator{T<:Number,B<:Operator} <: Operator{T}
    op::B
    uplo::Char
end

SymmetricOperator(B::Operator{T}, uplo::Symbol=:U) where {T<:Number} = SymmetricOperator{T,typeof(B)}(B, char_uplo(uplo))
convert(::Type{Operator{T}},A::SymmetricOperator) where {T}=SymmetricOperator(convert(Operator{T},A.op), A.uplo)

domainspace(P::SymmetricOperator)=domainspace(P.op)
rangespace(P::SymmetricOperator)=rangespace(P.op)
domain(P::SymmetricOperator)=domain(P.op)
bandwidths(P::SymmetricOperator) = (b = bandwidth(P.op, P.uplo == 'L' ? 1 : P.uplo == 'U' ? 2 : 0); (b, b))

function getindex(P::SymmetricOperator,k::Integer,j::Integer)
    if P.uplo == 'L'
        if j > k
            P.op[j,k]
        else
            P.op[k,j]
        end
    elseif P.uplo == 'U'
        if j < k
            P.op[j,k]
        else
            P.op[k,j]
        end
    end
end

Symmetric(A::Operator, uplo::Symbol)=SymmetricOperator(A, uplo)
Symmetric(A::Operator)=SymmetricOperator(A)


struct AdjointOperator{T<:Number,B<:Operator} <: Operator{T}
    op::B
end

AdjointOperator(B::Operator{T}) where {T<:Number}=AdjointOperator{T,typeof(B)}(B)
convert(::Type{Operator{T}},A::AdjointOperator) where {T}=AdjointOperator(convert(Operator{T},A.op))

domainspace(P::AdjointOperator)=rangespace(P.op)
rangespace(P::AdjointOperator)=domainspace(P.op)
domain(P::AdjointOperator)=domain(P.op)
bandwidths(P::AdjointOperator) = reverse(bandwidths(P.op))

getindex(P::AdjointOperator,k::Integer,j::Integer) = conj(P.op[j,k])

function BandedMatrix(S::SubOperator{T,TO}) where {T,TO<:AdjointOperator}
    kr,jr=parentindices(S)
    adjoint(BandedMatrix(view(parent(S).op,jr,kr)))
end

adjoint(A::Operator)=AdjointOperator(A)


struct TransposeOperator{T<:Number,B<:Operator} <: Operator{T}
    op::B
end

TransposeOperator(B::Operator{T}) where {T<:Number}=TransposeOperator{T,typeof(B)}(B)
convert(::Type{Operator{T}},A::TransposeOperator) where {T}=TransposeOperator(convert(Operator{T},A.op))

domainspace(P::TransposeOperator)=rangespace(P.op)
rangespace(P::TransposeOperator)=domainspace(P.op)
domain(P::TransposeOperator)=domain(P.op)
bandwidths(P::TransposeOperator) = reverse(bandwidths(P.op))

getindex(P::TransposeOperator,k::Integer,j::Integer) = P.op[j,k]

function BandedMatrix(S::SubOperator{T,TO}) where {T,TO<:TransposeOperator}
    kr,jr=parentindices(S)
    transpose(BandedMatrix(view(parent(S).op,jr,kr)))
end

transpose(A::Operator)=TransposeOperator(A)
