using ApproxFun
import Base:getindex

immutable RaggedMatrix{T} <: AbstractMatrix{T}
    data::Vector{T} # a Vector of non-zero entries
    cols::Vector{Int} # a Vector specifying the first index of each column
    m::Int #Number of rows
    function RaggedMatrix(data::Vector{T},cols::Vector{Int},m::Int)
        # make sure the cols are monitonically increasing
        @assert 1==cols[1]
        for j=1:length(cols)-1
            @assert cols[j] ≤ cols[j+1]
        end
        @assert cols[end] == length(data)+1

        # make sure we have less entries than the size of the matrix
        @assert length(data) ≤ m*(length(cols)-1)

        new(data,cols,m)
    end
end

RaggedMatrix(dat::Vector,cols::Vector{Int},m::Int) =
    RaggedMatrix{eltype(dat)}(dat,cols,m)

RaggedMatrix{T}(::Type{T},m::Int,colns::Vector{Int}) =
    RaggedMatrix(Array(T,sum(colns)),[1;cumsum(colns)],m)

RaggedMatrix(m::Int,collengths::Vector{Int}) = RaggedMatrix(Float64,m,collengths)

Base.size(A::RaggedMatrix) = (A.m,length(A.cols)-1)

Base.linearindexing{RM<:RaggedMatrix}(::Type{RM}) = Base.LinearSlow()

function getindex(A::RaggedMatrix,k::Int,j::Int)
    if k>size(A,1) || k < 1 || j>size(A,2) || j < 1
        throw(BoundsError(A,(k,j)))
    end

    if A.cols[j]+k-1 < A.cols[j+1]
        A.data[A.cols[j]+k-1]
    else
        zero(eltype(A))
    end
end

function Base.setindex!(A::RaggedMatrix,v,k::Int,j::Int)
    if k>size(A,1) || k < 1 || j>size(A,2) || j < 1
        throw(BoundsError(A,(k,j)))
    end

    if A.cols[j]+k-1 < A.cols[j+1]
        A.data[A.cols[j]+k-1]
    else
        throw(BoundsError(A,(k,j)))
    end
end


for (op,bop) in ((:(Base.rand),:rrand),(:(Base.zeros),:rzeros),(:(Base.ones),:rones))
    @eval begin
        $bop{T}(::Type{T},m::Int,colns::AbstractVector{Int}) =
            RaggedMatrix($op(T,sum(colns)),[1;1+cumsum(colns)],m)
        $bop(m::Int,colns::AbstractVector{Int}) = $bop(Float64,m,colns)
    end
end


colns=[1,2,3]
sum(colns),[1;cumsum(colns)]

rrand(10,1:10)
