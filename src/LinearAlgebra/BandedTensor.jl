
# Represents a BandedTensor.  In this case a Tensor is identified with
#  a block banded matrix.


abstract AbstractTensor{T} <: AbstractMatrix{T}

getindex(A::AbstractTensor,k::Integer,j::Integer) =



immutable BandedTensor{T} <: AbstractTensor{T}
    data::Vector{T}
    l::Int  # block lower bandwidth
    u::Int  # block upper bandwidth
    λ::Int  # sub lower bandwidth
    μ::Int  # sub upper bandwidth
end




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
