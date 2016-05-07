export FiniteOperator




immutable FiniteOperator{AT<:AbstractMatrix,T<:Number} <: BandedOperator{T}
    matrix::AT
end


FiniteOperator(M::AbstractMatrix)=FiniteOperator{typeof(M),eltype(M)}(M)


getindex(F::FiniteOperator,k::Integer,j::Integer) =
    k ≤ size(F.matrix,1) && j ≤ size(F.matrix,2) ? F.matrix[k,j] : zero(eltype(F))

function Base.copy{AT<:BandedMatrix,T}(S::SubBandedMatrix{T,FiniteOperator{AT,T}})
    kr,jr=parentindexes(S)
    if last(kr[1]) ≤ size(S.matrix,1) &&
        last(jr[2]) ≤ size(S.matrix,2)
        matrix[kr,jr]
    else
        default_copy(S)
    end
end



bandinds(T::FiniteOperator)=(1-size(T.matrix,1),size(T.matrix,2)-1)
bandinds{AT<:BandedMatrix}(T::FiniteOperator{AT})=bandinds(T.matrix)


Base.maximum(K::FiniteOperator)=maximum(K.matrix)


immutable FiniteFunctional{S,T} <: Functional{T}
    data::Vector{T}
    domainspace::S
end

Base.convert{T}(::Type{Functional{T}},S::FiniteFunctional)=FiniteFunctional(convert(Vector{T},S.data),S.domainspace)

domainspace(S::FiniteFunctional)=S.domainspace
datalength(S::FiniteFunctional)=length(S.data)


Base.getindex{S,T}(B::FiniteFunctional{S,T},k::Integer)=k≤datalength(B)?B.data[k]:zero(T)
Base.getindex{S,T}(B::FiniteFunctional{S,T},kr::Range)=T[B[k] for k=kr]
