immutable SubMatrix{T,B,I} <: AbstractMatrix{T}
    parent::B
    indexes::I
    dims::Tuple{Int,Int}
end

SubMatrix(A,inds,d) = SubMatrix{eltype(A),typeof(A),typeof(inds)}(A,inds,d)

immutable SubBandedMatrix{T,B,I} <: AbstractBandedMatrix{T}
    parent::B
    indexes::I
    dims::Tuple{Int,Int}
    l::Int
    u::Int
end

SubBandedMatrix(A,inds,d,l,u) =
    SubBandedMatrix{eltype(A),typeof(A),typeof(inds)}(A,inds,d,l,u)


bandwidth(V::SubBandedMatrix,k) = k==1?V.l:V.u


view(A::Operator,kr::Range,jr::Range) =
    SubMatrix(A,(kr,jr),(length(kr),length(jr)))

function view(A::Operator,kr::Range,jr::Range)
    st=step(kr)
    if isbanded(A) && st == step(jr)
        kr1=first(kr)
        jr1=first(jr)
        l,u=(bandinds(A,1)+kr1-jr1)÷st,(bandinds(A,2)+kr1-jr1)÷st
        SubBandedMatrix(A,(kr,jr),(length(kr),length(jr)),-l,u)
    else
        SubMatrix(A,(kr,jr),(length(kr),length(jr)))
    end
end


function view(A::Operator,kr::UnitRange,jr::UnitRange)
    @assert isbanded(A)
    shft=first(kr)-first(jr)
    l,u=max(bandwidth(A,1)-shft,0),max(bandinds(A,2)+shft,0)
    SubBandedMatrix(A,(kr,jr),(length(kr),length(jr)),l,u)
end

view(A::SubBandedMatrix,kr::UnitRange,jr::UnitRange) =
    view(A.parent,A.indexes[1][kr],A.indexes[2][jr])



## BLAS and matrix routines
# We assume that copy may be overriden

BLAS.axpy!(a,X::SubBandedMatrix,Y::AbstractMatrix) = BLAS.axpy!(a,copy(X),Y)

# this is for operators that implement copy via axpy!

copy_axpy!(S::SubBandedMatrix) =
    BLAS.axpy!(1.0,S,bzeros(eltype(S),size(S,1),size(S,2),bandwidth(S,1),bandwidth(S,2)))







#
# function defaultaxpy!(a,X::SubBandedMatrix,Y::BandedMatrix)
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
# function defaultaxpy!(a,X::SubBandedMatrix,Y)
#      @assert size(X)==size(Y)
#
#      for (k,j) in eachbandedindex(X)
#          Y[k,j]+=a*X[k,j]
#      end
#
#      Y
# end


## SubBandedOperator

immutable SubBandedOperator{T,B,I} <: Operator{T}
    parent::B
    indexes::I
    l::Int
    u::Int
end

SubBandedOperator(A,inds,l,u) = SubBandedOperator{eltype(A),typeof(A),typeof(inds)}(A,inds,l,u)

function view(A::Operator,kr::AbstractCount,jr::AbstractCount)
    @assert isbanded(A)
    st=step(kr)
    @assert st==step(jr)  # Otherwise, its not a banded operator
    kr1=first(kr)
    jr1=first(jr)
    l,u=(bandinds(A,1)+kr1-jr1)÷st,(bandinds(A,2)+kr1-jr1)÷st
    SubBandedOperator(A,(kr,jr),-l,u)
end

bandwidth(S::SubBandedOperator,k::Integer) = ifelse(k==1,S.l,S.u)
bandinds(S::SubBandedOperator) = (-S.l,S.u)




for TYP in (:SubBandedMatrix,:SubMatrix,:SubBandedOperator)
    @eval begin
        size(V::$TYP) = V.dims
        unsafe_getindex(V::$TYP,k::Integer,j::Integer) = V.parent[V.indexes[1][k],V.indexes[2][j]]
        getindex(V::$TYP,k::Integer,j::Integer) = V.parent[V.indexes[1][k],V.indexes[2][j]]
        getindex(V::$TYP,k::Integer,j::Range) = V.parent[V.indexes[1][k],V.indexes[2][j]]
        getindex(V::$TYP,k::Range,j::Integer) = V.parent[V.indexes[1][k],V.indexes[2][j]]
        getindex(V::$TYP,k::Range,j::Range) = V.parent[V.indexes[1][k],V.indexes[2][j]]
        Base.parent(S::$TYP) = S.parent
        Base.parentindexes(S::$TYP) = S.indexes
    end
end


function view(A::SubBandedOperator,kr::UnitRange,jr::UnitRange)
    view(A.parent,A.indexes[1][kr],A.indexes[2][jr])
end
