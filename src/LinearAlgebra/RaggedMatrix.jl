# FiniteRange gives the nonzero entries in a row/column
struct FiniteRange end

getindex(A::AbstractMatrix,::Type{FiniteRange},j::Integer) = A[1:colstop(A,j),j]
getindex(A::AbstractMatrix,k::Integer,::Type{FiniteRange}) = A[k,1:rowstop(A,k)]

const ⤓ = FiniteRange

mutable struct RaggedMatrix{T} <: AbstractMatrix{T}
    data::Vector{T} # a Vector of non-zero entries
    cols::Vector{Int} # a Vector specifying the first index of each column
    m::Int #Number of rows
    function RaggedMatrix{T}(data::Vector{T},cols::Vector{Int},m::Int) where T
        # make sure the cols are monitonically increasing
        @assert 1==cols[1]
        for j=1:length(cols)-1
            @assert cols[j] ≤ cols[j+1]
        end
        @assert cols[end] == length(data)+1

        # make sure we have less entries than the size of the matrix
        @assert length(data) ≤ m*(length(cols)-1)

        new{T}(data,cols,m)
    end
end

RaggedMatrix(dat::Vector,cols::Vector{Int},m::Int) =
    RaggedMatrix{eltype(dat)}(dat,cols,m)

RaggedMatrix(::Type{T},m::Int,colns::AbstractVector{Int}) where {T} =
    RaggedMatrix(Vector{T}(sum(colns)),Int[1;1+cumsum(colns)],m)

RaggedMatrix(m::Int,collengths::AbstractVector{Int}) = RaggedMatrix(Float64,m,collengths)


Base.size(A::RaggedMatrix) = (A.m,length(A.cols)-1)

colstart(A::RaggedMatrix,j::Integer) = 1
colstop(A::RaggedMatrix,j::Integer) = min(A.cols[j+1]-A.cols[j],size(A,1))

Base.IndexStyle(::Type{RM}) where {RM<:RaggedMatrix} = IndexCartesian()

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
        A.data[A.cols[j]+k-1]=v
    elseif v ≠ 0
        throw(BoundsError(A,(k,j)))
    end
    v
end

convert(::Type{RaggedMatrix{T}},R::RaggedMatrix{T}) where T = R

convert(::Type{RaggedMatrix{T}},R::RaggedMatrix) where T =
    RaggedMatrix{T}(Vector{T}(R.data), R.cols, R.m)


function convert(::Type{Matrix{T}},A::RaggedMatrix) where T
    ret = zeros(T,size(A,1),size(A,2))
    for j=1:size(A,2)
        ret[1:colstop(A,j),j] = view(A,1:colstop(A,j),j)
    end
    ret
end

convert(::Type{Matrix},A::RaggedMatrix) = Matrix{eltype(A)}(A)

Base.full(A::RaggedMatrix) = convert(Matrix,A)

function convert(::Type{RaggedMatrix{T}},B::BandedMatrix) where T
    l = bandwidth(B,1)
    ret = rzeros(T,size(B,1),Int[colstop(B,j) for j=1:size(B,2)])
    for j=1:size(B,2),k=colrange(B,j)
        ret[k,j] = B[k,j]
    end
    ret
end

convert(::Type{RaggedMatrix},B::BandedMatrix) = RaggedMatrix{eltype(B)}(B)

function convert(::Type{RaggedMatrix{T}},B::AbstractMatrix) where T
    ret = rzeros(T,size(B,1),Int[colstop(B,j) for j=1:size(B,2)])
    for j=1:size(B,2),k=colrange(B,j)
        ret[k,j] = B[k,j]
    end
    ret
end

convert(::Type{RaggedMatrix},B::AbstractMatrix) = RaggedMatrix{eltype(B)}(B)

Base.similar(B::RaggedMatrix,::Type{T}) where {T} = RaggedMatrix(Vector{T}(length(B.data)),copy(B.cols),B.m)

for (op,bop) in ((:(Base.rand),:rrand),(:(Base.zeros),:rzeros),(:(Base.ones),:rones))
    @eval begin
        $bop(::Type{T},m::Int,colns::AbstractVector{Int}) where {T} =
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
        for j = 1:size(X,2)
            @assert colstop(X,j) ≤ colstop(Y,j)  # check zeros otherwise
        end

        for j = 1:size(X,2)
            cs = colstop(X,j)
            BLAS.axpy!(a,view(X.data,X.cols[j]:X.cols[j]+cs-1),
                         view(Y.data,Y.cols[j]:Y.cols[j]+cs-1))
        end
    end
    Y
end

colstop(X::SubArray{T,2,RaggedMatrix{T},Tuple{UnitRange{Int},UnitRange{Int}}},
     j::Integer) where {T} = min(colstop(parent(X),j + first(parentindexes(X)[2])-1) -
                                            first(parentindexes(X)[1]) + 1,
                            size(X,1))

function BLAS.axpy!(a,X::RaggedMatrix,
                    Y::SubArray{T,2,RaggedMatrix{T},Tuple{UnitRange{Int},UnitRange{Int}}}) where T
    if size(X) ≠ size(Y)
        throw(BoundsError())
    end

    P = parent(Y)
    ksh = first(parentindexes(Y)[1]) - 1  # how much to shift
    jsh = first(parentindexes(Y)[2]) - 1  # how much to shift

    for j=1:size(X,2)
        cx=colstop(X,j)
        cy=colstop(Y,j)
        if cx > cy
            for k=cy+1:cx
                if X[k,j] ≠ 0
                    throw(BoundsError("Trying to add a non-zero to a zero."))
                end
            end
            kr = X.cols[j]:X.cols[j]+cy-1
        else
            kr = X.cols[j]:X.cols[j+1]-1
        end


        BLAS.axpy!(a,view(X.data,kr),
                    view(P.data,(P.cols[j + jsh] + ksh-1) + (1:length(kr))))
    end

    Y
end


function *(A::RaggedMatrix,B::RaggedMatrix)
    cols = zeros(Int,size(B,2))
    T = promote_type(eltype(A),eltype(B))
    for j=1:size(B,2),k=1:colstop(B,j)
        cols[j] = max(cols[j],colstop(A,k))
    end

    unsafe_A_mul_B!(RaggedMatrix(T,size(A,1),cols),A,B)
end

function unsafe_A_mul_B!(Y::RaggedMatrix,A::RaggedMatrix,B::RaggedMatrix)
    fill!(Y.data,0)

    for j=1:size(B,2),k=1:colstop(B,j)
        BLAS.axpy!(B[k,j],view(A,1:colstop(A,k),k),view(Y.data,Y.cols[j]-1+(1:colstop(A,k))))
    end

    Y
end

function Base.A_mul_B!(Y::RaggedMatrix,A::RaggedMatrix,B::RaggedMatrix)
    for j=1:size(B,2)
        col = 0
        for k=1:colstop(B,j)
            col = max(col,colstop(A,k))
        end

        if col > colstop(Y,j)
            throw(BoundsError())
        end
    end

    unsafe_A_mul_B!(Y,A,B)
end
