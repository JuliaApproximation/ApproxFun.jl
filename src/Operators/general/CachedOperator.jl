

export cache





## CachedOperator

mutable struct CachedOperator{T<:Number,DM<:AbstractMatrix,M<:Operator,DS,RS,BI} <: Operator{T}
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
        CachedOperator(BandedMatrix,op;padding=padding)
    elseif isbandedblockbanded(op) && !padding
        CachedOperator(BandedBlockBandedMatrix,op)
    elseif isblockbanded(op)
        CachedOperator(BlockBandedMatrix,op;padding=padding)
    elseif israggedbelow(op)
        CachedOperator(RaggedMatrix,op;padding=padding)
    else
        CachedOperator(Matrix,op;padding=padding)
    end
end

CachedOperator(op::Operator;padding::Bool=false) = default_CachedOperator(op;padding=padding)

"""
    cache(op::Operator)

Caches the entries of an operator, to speed up multiplying a Fun by the operator.
"""
cache(O::Operator;kwds...) = CachedOperator(O;kwds...)
cache(::Type{MT},O::Operator;kwds...) where {MT<:AbstractMatrix} = CachedOperator(MT,O;kwds...)

convert(::Type{Operator{T}},S::CachedOperator{T}) where {T} = S
convert(::Type{Operator{T}},S::CachedOperator) where {T} =
    CachedOperator(convert(Operator{T}, S.op),convert(AbstractMatrix{T}, S.data),
                    S.datasize,S.domainspace,S.rangespace,S.bandinds)


domainspace(C::CachedOperator) = C.domainspace
rangespace(C::CachedOperator) = C.rangespace
bandinds(C::CachedOperator{T,BM,M}) where {T<:Number,BM<:BandedMatrix,M<:Operator} = C.bandinds
blockbandinds(C::CachedOperator{T,BM,M}) where {T<:Number,BM<:BandedMatrix,M<:Operator} = C.bandinds


# TODO: cache this information as well
@wrapperstructure CachedOperator


@propagate_inbounds function Base.getindex(B::CachedOperator,k::Integer,j::Integer)
    resizedata!(B,k,j)
    if k ≤ size(B.data,1) && j ≤ size(B.data,2)
        B.data[k,j]
    else
        zero(eltype(B))
    end
end

function Base.getindex(B::CachedOperator,k::AbstractRange,j::AbstractRange)
    if !isempty(k) && !isempty(j)
        resizedata!(B,maximum(k),maximum(j))
        B.data[k,j]
    else
        Matrix{eltype(B)}(length(k),length(j))
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

@propagate_inbounds function Base.setindex!(B::CachedOperator,v,k,j)
    resizedata!(B,k,j)
    B.data[k,j] = v
end

@propagate_inbounds function Base.setindex!(B::CachedOperator,v,k::Integer)
    if size(B,1)==1
        B[1,k] = v
    elseif size(B,2)==1
        B[k,1] = v
    else
        error("Not implemented")
    end
 end






resizedata!(B::CachedOperator,::Colon,m::Integer) = resizedata!(B,size(B,1),m)
resizedata!(B::CachedOperator,n::Integer,::Colon) = resizedata!(B,n,size(B,2))


function mul_coefficients(B::CachedOperator,v::AbstractVector{T}) where T<:Number
    resizedata!(B,:,length(v))

    B.data*pad(v,size(B.data,2))
end
