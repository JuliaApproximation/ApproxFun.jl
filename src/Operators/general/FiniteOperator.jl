export FiniteOperator




immutable FiniteOperator{AT<:AbstractMatrix,T<:Number,DS,RS} <: Operator{T}
    matrix::AT
    domainspace::DS
    rangespace::RS
end


FiniteOperator(M::AbstractMatrix,ds::Space,rs::Space) =
    FiniteOperator{typeof(M),eltype(M),typeof(ds),typeof(rs)}(M,ds,rs)

FiniteOperator(M::AbstractMatrix) = FiniteOperator(M,ℓ⁰,ℓ⁰)

Base.convert{T}(::Type{Operator{T}},F::FiniteOperator) =
    FiniteOperator(convert(AbstractMatrix{T},F.matrix),F.domainspace,F.rangespace)::Operator{T}

domainspace(F::FiniteOperator) = F.domainspace
rangespace(F::FiniteOperator) = F.rangespace

function getindex(F::FiniteOperator,k::Integer,j::Integer)
    if k ≤ size(F.matrix,1) && j ≤ size(F.matrix,2)
        F.matrix[k,j]
    elseif k ≤ size(F,1) && j ≤ size(F,2)
        zero(eltype(F))
    else
        throw(BoundsError())
    end
end


function getindex(F::FiniteOperator,k::Integer)
    @assert size(F,1) == 1
    if k ≤ length(F.matrix)
        F.matrix[k]
    else
        zero(eltype(F))
    end
end

function Base.convert{AT<:BandedMatrix,T}(::Type{BandedMatrix},S::SubOperator{T,FiniteOperator{AT,T}})
    kr,jr=parentindexes(S)
    if last(kr[1]) ≤ size(S.matrix,1) &&
        last(jr[2]) ≤ size(S.matrix,2)
        matrix[kr,jr]
    else
        default_copy(S)
    end
end



bandinds(T::FiniteOperator) = (1-size(T.matrix,1),size(T.matrix,2)-1)
bandinds{AT<:BandedMatrix}(T::FiniteOperator{AT}) = bandinds(T.matrix)


Base.maximum(K::FiniteOperator) = maximum(K.matrix)
