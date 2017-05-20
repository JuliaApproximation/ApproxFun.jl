export FiniteOperator




immutable FiniteOperator{AT<:AbstractMatrix,T<:Number,DS,RS} <: Operator{T}
    matrix::AT
    domainspace::DS
    rangespace::RS
end


FiniteOperator(M::AbstractMatrix,ds::Space,rs::Space) =
    FiniteOperator{typeof(M),eltype(M),typeof(ds),typeof(rs)}(M,ds,rs)

FiniteOperator(M::AbstractMatrix) =
    FiniteOperator(M,EuclideanSpace(size(M,2)),EuclideanSpace(size(M,1)))

Base.convert{T}(::Type{Operator{T}},F::FiniteOperator) =
    FiniteOperator(convert(AbstractMatrix{T},F.matrix),F.domainspace,F.rangespace)::Operator{T}


Base.promote_rule{OT<:Operator,MT<:AbstractMatrix}(::Type{OT},::Type{MT}) = Operator{promote_type(eltype(OT),eltype(MT))}

Base.convert{T}(::Type{Operator{T}},M::AbstractMatrix) = FiniteOperator(AbstractMatrix{T}(M))
Base.convert(::Type{Operator},M::AbstractMatrix) = Operator{eltype(M)}(M)

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


bandinds(T::FiniteOperator) = bandinds(T.matrix)
# TODO: What should the blockbandinds be?
Base.maximum(K::FiniteOperator) = maximum(K.matrix)
