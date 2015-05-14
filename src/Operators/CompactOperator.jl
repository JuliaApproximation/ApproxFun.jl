export CompactOperator



immutable CompactOperator{T<:Number} <: BandedOperator{T}
    matrix::Array{T,2}
end



function matrix_addentries!(M::Array,A,kr::Range)
    for k=kr[1]:min(size(M,1),kr[end]),j=1:size(M,2)
        A[k,j] += M[k,j]
    end
    
    A
end


addentries!(T::CompactOperator,A,kr::Range)=matrix_addentries!(T.matrix,A,kr)

bandinds(T::CompactOperator)=(1-size(T.matrix,1),size(T.matrix,2)-1)


immutable CompactFunctional{S,T} <: Functional{T}
    data::Vector{T}
    domainspace::S
end

Base.convert{T}(::Type{Functional{T}},S::CompactFunctional)=CompactFunctional(convert(Vector{T},S.data),S.domainspace)

domainspace(S::CompactFunctional)=S.domainspace
datalength(S::CompactFunctional)=length(S.data)


Base.getindex{S,T}(B::CompactFunctional{S,T},k::Integer)=k≤datalength(B)?B.data[k]:zero(T)
Base.getindex{S,T}(B::CompactFunctional{S,T},kr::Range)=T[B[k] for k=kr]
