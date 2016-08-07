

export cache





## CachedOperator

type CachedOperator{T<:Number,DM<:AbstractMatrix,M<:Operator,DS,RS,BI} <: Operator{T}
    op::M
    data::DM
    datasize::Tuple{Int,Int}
    domainspace::DS
    rangespace::RS
    bandinds::BI
end

CachedOperator(op::Operator,data::AbstractMatrix,sz::Tuple{Int,Int},ds,rs,bi) =
    CachedOperator{eltype(data),typeof(data),typeof(op),
                    typeof(ds),typeof(rs),typeof(bi)}(op,data,sz,ds,rs,bi)


CachedOperator(op::Operator,data::AbstractMatrix,sz::Tuple{Int,Int}) =
    CachedOperator(op,data,sz,domainspace(op),rangespace(op),bandinds(op))
CachedOperator(op::Operator,data::AbstractMatrix) = CachedOperator(op,data,size(data))
function CachedOperator(op::Operator)
    if isbanded(op)
        CachedOperator(op,BandedMatrix(eltype(op),0,0,bandwidth(op,1),bandwidth(op,2)))
    else
        CachedOperator(op,Array(eltype(op),0,0))
    end
end


cache(O::Operator) = CachedOperator(O)

Base.convert{T}(::Type{Operator{T}},S::CachedOperator{T}) = S
Base.convert{T}(::Type{Operator{T}},S::CachedOperator) =
    CachedOperator(Operator{T}(S.op),AbstractMatrix{T}(S.data),
                    S.datasize,S.domainspace,S.rangespace,S.bandinds)


domainspace(C::CachedOperator) = C.domainspace
rangespace(C::CachedOperator) = C.rangespace
bandinds(C::CachedOperator) = C.bandinds
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


function resizedata!{T<:Number}(B::CachedOperator{T,Matrix{T}},n::Integer,m::Integer)
    if n > size(B,1) || m > size(B,2)
        throw(ArgumentError("Cannot resize beyound size of operator"))
    end

    # this does nothing if already in dimensions
    N,M=size(B.data)
    if n > N && m > M
        B.data = unsafe_resize!(B.data,n,m)
    elseif n > N
        B.data = unsafe_resize!(B.data,n,:)
    elseif m > M
        B.data = unsafe_resize!(B.data,:,m)
    end

    if n ≤ B.datasize[1] && m ≤ B.datasize[2]
        # do nothing
        B
    elseif n ≤ B.datasize[1]
        kr,jr=1:B.datasize[1],B.datasize[2]+1:m
        B.data[kr,jr] = B.op[kr,jr]
        B.datasize = (B.datasize[1],m)
        B
    elseif m ≤ B.datasize[2]
        kr,jr=B.datasize[1]+1:n,1:B.datasize[2]
        B.data[kr,jr] = B.op[kr,jr]
        B.datasize = (n,B.datasize[2])
        B
    else
        # resize rows then columns
        resizedata!(resizedata!(B,n,B.datasize[2]),n,m)
    end
end

function resizedata!{T<:Number}(B::CachedOperator{T,BandedMatrix{T}},n::Integer,::Colon)
    if n > size(B,1)
        throw(ArgumentError("Cannot resize beyound size of operator"))
    end

    if n > B.datasize[1]
        pad!(B.data,2n,:)

        kr=B.datasize[1]+1:n
        jr=max(B.datasize[1]+1-B.data.l,1):n+B.data.u
        BLAS.axpy!(1.0,view(B.op,kr,jr),view(B.data,kr,jr))

        B.datasize = (n,n+B.data.u)
    end

    B
end


resizedata!{T<:Number}(B::CachedOperator{T,BandedMatrix{T}},n::Integer,m::Integer) =
    resizedata!(B,n,:)

resizedata!(B::CachedOperator,::Colon,m::Integer) = resizedata!(B,B.datasize[1],m)
resizedata!(B::CachedOperator,n::Integer,::Colon) = resizedata!(B,n,B.datasize[2])
