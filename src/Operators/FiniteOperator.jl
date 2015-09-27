export FiniteOperator




immutable FiniteOperator{T<:Number} <: BandedOperator{T}
    matrix::Array{T,2}
end



function matrix_addentries!(M::Array,A,kr::Range)
    for k=kr[1]:min(size(M,1),kr[end]),j=1:size(M,2)
        A[k,j] += M[k,j]
    end

    A
end


addentries!(T::FiniteOperator,A,kr::Range,::Colon)=matrix_addentries!(T.matrix,A,kr)

bandinds(T::FiniteOperator)=(1-size(T.matrix,1),size(T.matrix,2)-1)


# An infinite slice of a FiniteOperator is also a FiniteOperator

function Base.slice(K::FiniteOperator,kr::FloatRange,jr::FloatRange)
    st=step(kr)
    @assert step(jr)==st
    @assert last(kr)==last(jr)==Inf
    FiniteOperator(K.matrix[first(kr):st:end,first(jr):st:end])
end


immutable FiniteFunctional{S,T} <: Functional{T}
    data::Vector{T}
    domainspace::S
end

Base.convert{T}(::Type{Functional{T}},S::FiniteFunctional)=FiniteFunctional(convert(Vector{T},S.data),S.domainspace)

domainspace(S::FiniteFunctional)=S.domainspace
datalength(S::FiniteFunctional)=length(S.data)


Base.getindex{S,T}(B::FiniteFunctional{S,T},k::Integer)=k≤datalength(B)?B.data[k]:zero(T)
Base.getindex{S,T}(B::FiniteFunctional{S,T},kr::Range)=T[B[k] for k=kr]
