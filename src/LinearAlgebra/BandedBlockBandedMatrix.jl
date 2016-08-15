using BandedMatrices
import Base.BLAS.BlasFloat

# Represents a block banded matrix with banded blocks


function blocklookup(rows)
    rowblocks=Array(Int,0)
    for K in rows, k in 1:K
        push!(rowblocks,K)
    end
    rowblocks
end



immutable BandedBlockBandedMatrix{T,RI,CI} <: AbstractMatrix{T}
    data::Matrix{T}

    l::Int  # block lower bandwidth
    u::Int  # block upper bandwidth
    λ::Int  # sub lower bandwidth
    μ::Int  # sub upper bandwidth

    rows::RI   # the size of the row blocks
    cols::CI  # the size of the column blocks

    rowblocks::Vector{Int}
    colblocks::Vector{Int}

    function BandedBlockBandedMatrix(data,l,u,λ,μ,rows,cols)
        @assert size(data,1) == λ+μ+1
        new(data,l,u,λ,μ,rows,cols,blocklookup(rows),blocklookup(cols))
    end
end

BandedBlockBandedMatrix(data,l,u,λ,μ,rows,cols) =
    BandedBlockBandedMatrix{eltype(data),typeof(rows),typeof(cols)}(data,l,u,λ,μ,rows,cols)


Base.size(A::BandedBlockBandedMatrix) = sum(A.rows),sum(A.cols)
function getblock(A::BandedBlockBandedMatrix,K,J)
    if K < J-A.u || K > J+A.l
        zeros(eltype(A),A.rows[K],A.cols[J])
    else
        # column block K-J+A.u+1,J
        S=sum(A.cols[1:J-1])*(l+u+1)  # number of columns before current block
        cols=S+(K-J+A.u)A.cols[J]+1:S+(K-J+A.u+1)A.cols[J]
        BandedMatrix(A.data[:,cols],
                        A.rows[K],A.λ,A.μ)
    end
end

function viewblock{T<:BlasFloat}(A::BandedBlockBandedMatrix{T},K,J)
    if K < J-A.u || K > J+A.l
        error("Cannot view zero blocks")
    else
        # column block K-J+A.u+1,J
        S=sum(A.cols[1:J-1])*(l+u+1)  # number of columns before current block
        p=pointer(A.data)
        st=stride(A.data,2)
        sz=sizeof(p)
        cols=S+(K-J+A.u)A.cols[J]+1:S+(K-J+A.u+1)A.cols[J]

        p+=(first(cols)-1)*st*sz
        BandedMatrix(unsafe_wrap(Array,p,(st,length(cols))),
                        A.rows[K],A.λ,A.μ)
    end
end


function Base.getindex{T<:BlasFloat}(A::BandedBlockBandedMatrix{T},k,j)
    K=A.rowblocks[k];J=A.colblocks[j]
    if K < J-A.u || K > J+A.l
        zero(eltype(A))
    else
        k2=k-sum(rows[1:K-1])
        j2=j-sum(cols[1:J-1])
        viewblock(A,K,J)[k2,j2]
    end
end

function Base.setindex!(A::BandedBlockBandedMatrix,v,k,j)
    K=A.rowblocks[k];J=A.colblocks[j]
    k2=k-sum(rows[1:K-1])
    j2=j-sum(cols[1:J-1])
    setindex!(viewblock(A,K,J),v,k2,j2)
end

Base.linearindexing{BBBM<:BandedBlockBandedMatrix}(::Type{BBBM}) =
    Base.LinearSlow()



l=u=1
λ=μ=1
N=M=50
cols=rows=1:N
@time data=ones(λ+μ+1,(l+u+1)*sum(cols))
data=1.0reshape(1:length(data)|>collect,size(data,1),size(data,2))
@time A=BandedBlockBandedMatrix(data,l,u,λ,μ,rows,cols)
    rand(5,5)

A[1,1]

# skipped

colblocks=rowblocks

k,j=5,5


getblock(A,K,J)

K,J
rowblocks

K,J

k
j
k2,j2

getblock(A,K,J)

rowblocks
colblocks[j-1]

k,j=3,4
    sum(cols[1:colblocks[j-1]])

sum(1:2)
colblocks


sum(cols[1:colblocks[j-1]])

@time blocklookup(1:50)
    3


sum(1:500)

rows

rowblocks


getblock(A,3,1)


BandedMatrix(A.data[:,u*columns[1]+1:(u+1)columns[1]],A.rows[1],A.λ,A.μ)
BandedMatrix(A.data[:,(u+1)columns[1]+1:(u+2)columns[1]],A.rows[2],A.λ,A.μ)
BandedMatrix(A.data[:,(u+2)columns[1]+1:(u+3)columns[1]],A.rows[3],A.λ,A.μ)

(u-1)

K=3;J=2
BandedMatrix(A.data[:,(l+u+1)columns[1]+u*columns[2]+1:(l+u+1)columns[1]+(u+1)*columns[2]],
                    A.rows[2],A.λ,A.μ)


A.data

K-J+A.u

K-J+A.u+1

u-1

K=1;J=2



A.data


columns
sum(columns)

length(data)

sum(rows)^2


sum(rows)

reshape(collect(1:16),2,2,2,2)


KO=C.op

kr=1:5
jr=1:5
N=K=maximum(kr)
M=J=maximum(jr)
A=KO.ops[1][1:K,1:J]
B=KO.ops[2][1:N,1:M]

KR=kron(full(B),full(A))
findfirst(KO.domaintensorizer,(1,4))
KR[1,findfirst(KO.domaintensorizer,(2,1))]

A[1,1]*B[1,3]
KO[1:25,1:25]

kron2tensor(A
