

export cache





## CachedOperator

type CachedOperator{T<:Number,DM<:AbstractMatrix,M<:Operator,DS,RS,BI} <: Operator{T}
    op::M
    data::DM
    datasize::Tuple{Int,Int}
    domainspace::DS
    rangespace::RS
    bandinds::BI
    padding::Bool   # records whether the operator should be padded for QR operations
                    #  note that BandedMatrix/AlmostBandedMatrix don't need this
                    #  as the padding is done in the bandinds
end

CachedOperator(op::Operator,data::AbstractMatrix,sz::Tuple{Int,Int},ds,rs,bi,pd=false) =
    CachedOperator{eltype(data),typeof(data),typeof(op),
                    typeof(ds),typeof(rs),typeof(bi)}(op,data,sz,ds,rs,bi,pd)


CachedOperator(op::Operator,data::AbstractMatrix,sz::Tuple{Int,Int},pd=false) =
    CachedOperator(op,data,sz,domainspace(op),rangespace(op),bandinds(op),pd)
CachedOperator(op::Operator,data::AbstractMatrix,padding=false) = CachedOperator(op,data,size(data),padding)


function default_CachedOperator(op::Operator;padding::Bool=false)
    if isbanded(op)
        l,u=bandwidths(op)
        padding && (u+=l)
        data=BandedMatrix(eltype(op),0,0,l,u)
        CachedOperator(op,data,size(data),domainspace(op),rangespace(op),(-l,u),padding)
    elseif isbandedblock(op)
        l,u=blockbandwidths(op)
        padding=false
        padding && (u+=l)
        data=BandedBlockMatrix(eltype(op),l,u,1:0,1:0)  # TODO: type of rows/cols
        CachedOperator(op,data,size(data),domainspace(op),rangespace(op),(-l,u),padding)
    elseif israggedbelow(op)
        CachedOperator(op,RaggedMatrix(eltype(op),0,Int[]),padding)
    else
        CachedOperator(op,Array(eltype(op),0,0),padding)
    end
end

CachedOperator(op::Operator;padding::Bool=false) = default_CachedOperator(op;padding=padding)


cache(O::Operator;kwds...) = CachedOperator(O;kwds...)

Base.convert{T}(::Type{Operator{T}},S::CachedOperator{T}) = S
Base.convert{T}(::Type{Operator{T}},S::CachedOperator) =
    CachedOperator(Operator{T}(S.op),AbstractMatrix{T}(S.data),
                    S.datasize,S.domainspace,S.rangespace,S.bandinds)


domainspace(C::CachedOperator) = C.domainspace
rangespace(C::CachedOperator) = C.rangespace
bandinds{T,BM<:BandedMatrix}(C::CachedOperator{T,BM}) = C.bandinds
blockbandinds{T,BM<:BandedBlockMatrix}(C::CachedOperator{T,BM}) = C.bandinds
Base.stride(C::CachedOperator) = stride(C.op)

# when bandinds are infinite, use colstart/colstop from operator
# TODO: cache this information as well
for func in (:(ApproxFun.colstart),:(ApproxFun.colstop),
                :(ApproxFun.rowstart),:(ApproxFun.rowstop))
    @eval $func{T,DM,M,DS,RS}(D::CachedOperator{T,DM,M,DS,RS,Tuple{Infinity{Bool},Infinity{Bool}}},
                         k::Integer) = $func(D.op,k)
end


function Base.getindex(B::CachedOperator,k::Integer,j::Integer)
    resizedata!(B,k,j)
    B.data[k,j]
end

function Base.getindex(B::CachedOperator,k::Range,j::Range)
    if !isempty(k) && !isempty(j)
        resizedata!(B,maximum(k),maximum(j))
        B.data[k,j]
    else
        Array(eltype(B),length(k),length(j))
    end
end

function Base.getindex(B::CachedOperator,k::Integer)
    if size(B,1)==1
        B[1,k]
    elseif size(B,2)==1
        B[k,1]
    else
        error("Not implemented")
    end
 end


resizedata!(B::CachedOperator,::Colon,m::Integer) = resizedata!(B,size(B,1),m)
resizedata!(B::CachedOperator,n::Integer,::Colon) = resizedata!(B,n,size(B,2))
