

export SavedBandedOperator, cache





## CachedOperator

type CachedOperator{T<:Number,M<:Operator,DS,RS,BI} <: Operator{T}
    op::M
    data::Matrix{T}
    datasize::Tuple{Int,Int}
    domainspace::DS
    rangespace::RS
    bandinds::BI
end

CachedOperator(op::Operator,data,sz) =
    CachedOperator(op,data,sz,domainspace(op),rangespace(op),bandinds(op))
CachedOperator(op::Operator,data) = CachedOperator(op,data,size(data))
CachedOperator(op::Operator) = CachedOperator(op,Array(eltype(op),0,0))


cache(O::Operator) = CachedOperator(O)

Base.convert{T}(::Type{Operator{T}},S::CachedOperator) =
    CachedOperator(Operator{T}(S.op),Matrix{T}(S.data),
                    S.datasize,S.domainspace,S.rangespace,S.bandinds)


domainspace(C::CachedOperator) = C.domainspace
rangespace(C::CachedOperator) = C.rangespace
bandinds(C::CachedOperator) = C.bandinds
Base.stride(C::CachedOperator) = stride(C.op)

# when bandinds are infinite, use colstart/colstop from operator
# TODO: cache this information as well
for func in (:(ApproxFun.colstart),:(ApproxFun.colstop),
                :(ApproxFun.rowstart),:(ApproxFun.rowstop))
    @eval $func{T,M,DS,RS}(D::CachedOperator{T,M,DS,RS,Tuple{Infinity{Bool},Infinity{Bool}}},
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


function resizedata!(B::CachedOperator,n::Integer,m::Integer)
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
        kr,jr=1:B.datasize[1]+1:n,1:B.datasize[2]
        B.data[kr,jr] = B.op[kr,jr]
        B.datasize = (n,B.datasize[2])
        B
    else
        # resize rows then columns
        resizedata!(resizedata!(B,n,B.datasize[2]),n,m)
    end
end

resizedata!(B::CachedOperator,::Colon,m::Integer) = resizedata!(B,B.datasize[1],m)
resizedata!(B::CachedOperator,n::Integer,::Colon) = resizedata!(B,n,B.datasize[2])



## SavedBandedOperator


type SavedBandedOperator{T<:Number,M<:Operator} <: Operator{T}
    op::M
    data::BandedMatrix{T}   #Shifted to encapsolate bandedness
    datalength::Int
    bandinds::Tuple{Int,Int}
end



# convert needs to throw away calculated data
function Base.convert{T}(::Type{Operator{T}},S::SavedBandedOperator)
    if eltype(S) == T
        S
    else
        SavedBandedOperator(convert(Operator{T},S.op),
                            convert(BandedMatrix{T},S.data),
                            S.datalength,S.bandinds)
    end
end


#TODO: index(op) + 1 -> length(bc) + index(op)
function SavedBandedOperator{T<:Number}(op::Operator{T})
    data = bzeros(T,0,:,bandinds(op))  # bzeros is needed to allocate top of array
    SavedBandedOperator(op,data,0,bandinds(op))
end



for OP in (:domain,:domainspace,:rangespace,:(Base.stride))
    @eval $OP(S::SavedBandedOperator)=$OP(S.op)
end

bandinds(B::SavedBandedOperator)=B.bandinds




function Base.getindex(B::SavedBandedOperator,k::Integer,j::Integer)
    resizedata!(B,k)
    B.data[k,j]
end

function resizedata!(B::SavedBandedOperator,n::Integer)
    if n > B.datalength
        pad!(B.data,2n,:)

        kr=B.datalength+1:n
        jr=max(B.datalength+1-B.data.l,1):n+B.data.u
        BLAS.axpy!(1.0,view(B.op,kr,jr),view(B.data,kr,jr))

        B.datalength = n
    end

    B
end
