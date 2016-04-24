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
        unsafe_getindex(V::$TYP,k::Integer,j::Integer) = V.parent[V.indexes[1][k],V.indexes[2][j]]
        getindex(V::$TYP,k::Integer,j::Integer) = V.parent[V.indexes[1][k],V.indexes[2][j]]
        getindex(V::$TYP,k,j) = V.parent[V.indexes[1][k],V.indexes[2][j]]
        Base.parent(S::$TYP) = S.parent
        Base.parentindexes(S::$TYP) = S.indexes
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



## BLAS and matrix routines
# We assume that copy may be overriden

BLAS.axpy!(a,X::SubBandedOperator,Y::AbstractMatrix) = axpy!(a,copy(X),Y)

# this is for operators that implement copy via axpy!

copy_axpy!(S::SubBandedOperator) =
    BLAS.axpy!(1.0,S,bzeros(eltype(S),size(S,1),size(S,2),bandwidth(S,1),bandwidth(S,2)))







#
# function defaultaxpy!(a,X::SubBandedOperator,Y::BandedMatrix)
#      @assert size(X)==size(Y)
#      @assert bandwidth(X,1) ≤ bandwidth(Y,1) && bandwidth(X,2) ≤ bandwidth(Y,2)
#
#      for (k,j) in eachbandedindex(X)
#          Y[k,j]+=a*X[k,j]
#      end
#
#      Y
# end
#
# function defaultaxpy!(a,X::SubBandedOperator,Y)
#      @assert size(X)==size(Y)
#
#      for (k,j) in eachbandedindex(X)
#          Y[k,j]+=a*X[k,j]
#      end
#
#      Y
# end
