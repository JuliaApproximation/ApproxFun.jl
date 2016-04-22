immutable SubOperator{T,B,I} <: AbstractMatrix{T}
    parent::B
    indexes::I
    dims::Tuple{Int,Int}
end


immutable SubBandedOperator{T,B,I} <: AbstractBandedMatrix{T}
    parent::B
    indexes::I
    dims::Tuple{Int,Int}
    l::Int
    u::Int
end

for TYP in (:SubBandedOperator,:SubOperator)
    @eval begin
        size(V::$TYP) = V.dims
        unsafe_getindex(V::$TYP,k::Integer,j::Integer)=V.parent[V.indexes[1][k],V.indexes[2][j]]
        getindex(V::$TYP,k::Integer,j::Integer)=V.parent[V.indexes[1][k],V.indexes[2][j]]
        getindex(V::$TYP,k,j)=V.parent[V.indexes[1][k],V.indexes[2][j]]
    end
end

bandwidth(V::SubBandedOperator,k)=k==1?V.l:V.u

Base.sub(A::Operator,kr::Range,jr::Range)=SubOperator{eltype(A),typeof(A),Tuple{typeof(kr),typeof(jr)}}(A,(kr,jr),(length(kr),length(jr)))


function Base.sub(A::BandedOperator,kr::UnitRange,jr::UnitRange)
    shft=first(kr)-first(jr)
    SubBandedOperator{eltype(A),typeof(A),Tuple{typeof(kr),typeof(jr)}}(A,(kr,jr),
                                                                      (length(kr),length(jr)),
                                                                      max(bandwidth(A,1)-shft,0),max(bandinds(A,2)+shft,0))
end
