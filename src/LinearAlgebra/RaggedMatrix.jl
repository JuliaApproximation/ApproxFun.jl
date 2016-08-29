type RaggedMatrix{T} <: AbstractMatrix{T}
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

colstop(A::RaggedMatrix,j::Integer) = A.cols[j+1]-A.cols[j]

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


## BLAS

function Base.A_mul_B!(y::Vector,A::RaggedMatrix,b::Vector)
    m=size(A,2)

    if m ≠ length(b) || size(A,1) ≠ length(y)
        throw(BoundsError())
    end
    T=eltype(y)
    fill!(y,zero(T))
    for j=1:m
        kr=A.cols[j]:A.cols[j+1]-1
        BLAS.axpy!(b[j],view(A.data,kr),view(y,1:length(kr)))
    end
    y
end


function BLAS.axpy!(a,X::RaggedMatrix,Y::RaggedMatrix)
    if size(X) ≠ size(Y)
        throw(BoundsError())
    end

    if X.cols == Y.cols
        BLAS.axpy!(a,X.data,Y.data)
    else
        error("Not implemented.")
    end
end

colstop{T}(X::SubArray{T,2,RaggedMatrix{T},Tuple{UnitRange{Int},UnitRange{Int}}},
        j::Integer) = min(colstop(parent(X),j + first(parentindexes(X)[2])-1) -
                                            first(parentindexes(X)[1]) + 1,
                            size(X,1))

function BLAS.axpy!{T}(a,X::RaggedMatrix,
                       Y::SubArray{T,2,RaggedMatrix{T},Tuple{UnitRange{Int},UnitRange{Int}}})
    if size(X) ≠ size(Y)
        throw(BoundsError())
    end

    for j=1:size(X,2)
        @assert colstop(X,j) == colstop(Y,j)
    end

    P = parent(Y)
    ksh = first(parentindexes(Y)[1]) - 1  # how much to shift
    jsh = first(parentindexes(Y)[2]) - 1  # how much to shift

    for j=1:size(X,2)
        kr = X.cols[j]:X.cols[j+1]-1
        BLAS.axpy!(a,view(X.data,kr),
                    view(P.data,(P.cols[j + jsh] + ksh-1) + (1:length(kr))))
    end

    Y
end
