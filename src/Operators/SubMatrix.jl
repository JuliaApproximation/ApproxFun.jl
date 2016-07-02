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

function view(A::BandedOperator,kr::Range,jr::Range)
    st=step(kr)
    if st == step(jr)
        kr1=first(kr)
        jr1=first(jr)
        l,u=(bandinds(A,1)+kr1-jr1)÷st,(bandinds(A,2)+kr1-jr1)÷st
        SubBandedMatrix(A,(kr,jr),(length(kr),length(jr)),-l,u)
    else
        SubMatrix(A,(kr,jr),(length(kr),length(jr)))
    end
end


function view(A::BandedOperator,kr::UnitRange,jr::UnitRange)
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

immutable SubBandedOperator{T,B,I} <: BandedOperator{T}
    parent::B
    indexes::I
    l::Int
    u::Int
end

SubBandedOperator(A,inds,l,u) = SubBandedOperator{eltype(A),typeof(A),typeof(inds)}(A,inds,l,u)

function view(A::BandedOperator,kr::AbstractCount,jr::AbstractCount)
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



#
# # Some of this is verbatim from IndexSlice
# immutable SliceOperator{T,B} <: BandedOperator{T}
#     op::B
#     rowindex::Int
#     colindex::Int
#     rowstride::Int
#     colstride::Int
#
#     function SliceOperator(o,r,c,rs,cs)
#         @assert rs == cs
#         @assert rs != 0
#         @assert mod(r-c,rs)==0
#         @assert mod(stride(o),rs)==0
#
#         new(o,r,c,rs,cs)
#     end
# end
#
# SliceOperator{T<:Number}(B::Operator{T},r,c,rs,cs)=SliceOperator{T,typeof(B)}(B,r,c,rs,cs)
# SliceOperator{T<:Number}(B::Operator{T},r,c,rs)=SliceOperator{T,typeof(B)}(B,r,c,rs,rs)
# SliceOperator{T<:Number}(B::Operator{T},r,c)=SliceOperator{T,typeof(B)}(B,r,c,1,1)
#
#
# Base.convert{BT<:Operator}(::Type{BT},S::SliceOperator)=SliceOperator(convert(BandedOperator{eltype(BT)},S.op),
#                                                                         S.rowindex,S.colindex,S.rowstride,S.colstride)
#
# bandinds(S::SliceOperator)=(div(bandinds(S.op,1)+S.rowindex-S.colindex,S.rowstride),div(bandinds(S.op,2)+S.rowindex-S.colindex,S.rowstride))
#
# function destride_addentries!(op,ri,ci,rs,cs,A,kr::UnitRange)
#     r1=rs*kr[1]+ri:rs:rs*kr[end]+ri
#     addentries!(op,IndexSlice(A,ri,ci,rs,cs),r1,:)
#     A
# end
#
# function destride_addentries!(op,ri,ci,A,kr::UnitRange)
#     r1=kr[1]+ri:kr[end]+ri
#     addentries!(op,IndexSlice(A,ri,ci,1,1),r1,:)
#     A
# end
#
# function destride_addentries!(S::SliceOperator,A,kr::Range)
#     if S.rowstride==S.colstride==1
#         destride_addentries!(S.op,S.rowindex,S.colindex,A,kr)
#     else
#         destride_addentries!(S.op,S.rowindex,S.colindex,S.rowstride,S.colstride,A,kr)
#     end
# end
#
# addentries!(S::SliceOperator,A,kr,::Colon)=destride_addentries!(S,A,kr)
# domain(S::SliceOperator)=domain(S.op)
# domainspace(S::SliceOperator)=S.colindex==0&&S.colstride==1?domainspace(S.op):SliceSpace(domainspace(S.op),S.colindex,S.colstride)
# rangespace(S::SliceOperator)=SliceSpace(rangespace(S.op),S.rowindex,S.rowstride)
