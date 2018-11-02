export FiniteOperator




struct FiniteOperator{AT<:AbstractMatrix,T<:Number,DS,RS} <: Operator{T}
    matrix::AT
    domainspace::DS
    rangespace::RS
end


FiniteOperator(M::AbstractMatrix{<:Number},ds::Space,rs::Space) =
    FiniteOperator{typeof(M),eltype(M),typeof(ds),typeof(rs)}(M,ds,rs)

FiniteOperator(M::AbstractMatrix{<:Number}) =
    FiniteOperator(M,EuclideanSpace(size(M,2)),EuclideanSpace(size(M,1)))

convert(::Type{Operator{T}},F::FiniteOperator) where {T} =
    FiniteOperator(convert(AbstractMatrix{T},F.matrix),F.domainspace,F.rangespace)::Operator{T}


Base.promote_rule(::Type{OT},::Type{MT}) where {OT<:Operator,MT<:AbstractMatrix} = Operator{promote_type(eltype(OT),eltype(MT))}

convert(::Type{Operator{T}},M::AbstractMatrix{<:Number}) where {T} = FiniteOperator(AbstractMatrix{T}(M))
convert(::Type{Operator},M::AbstractMatrix{<:Number}) = Operator{eltype(M)}(M)

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

function BandedMatrix(S::SubOperator{T,FiniteOperator{AT,T}}) where {AT<:BandedMatrix,T}
    kr,jr=parentindices(S)
    if last(kr[1]) ≤ size(S.matrix,1) &&
        last(jr[2]) ≤ size(S.matrix,2)
        matrix[kr,jr]
    else
        default_copy(S)
    end
end


bandwidths(T::FiniteOperator) = bandwidths(T.matrix)
blockbandwidths(T::FiniteOperator) = block(rangespace(T),size(T.matrix,1)).n[1]-1,block(domainspace(T),size(T.matrix,2)).n[1]-1
Base.maximum(K::FiniteOperator) = maximum(K.matrix)
