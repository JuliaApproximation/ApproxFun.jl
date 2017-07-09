# This file is based on rowvector.jl in Julia. License is MIT: https://julialang.org/license
# The motivation for this file is to allow RowVector which doesn't transpose the entries

import Base: convert, similar, length, size, indices, IndexStyle,
            IndexLinear, @propagate_inbounds, getindex, setindex!,
            broadcast, hcat, typed_hcat, A_mul_Bt, At_mul_Bt, At_mul_B,
            A_mul_Bc, Ac_mul_Bc, Ac_mul_B, At_ldiv_B, map, parent

import Base.LinAlg: check_types, check_tail_indices, to_vec

"""
    RowVector(vector)

A lazy-view wrapper of an `AbstractVector`, which turns a length-`n` vector into a `1×n`
shaped row vector and represents the transpose of a vector (the elements are also transposed
recursively). This type is usually constructed (and unwrapped) via the [`transpose`](@ref)
function or `.'` operator (or related [`ctranspose`](@ref) or `'` operator).

By convention, a vector can be multiplied by a matrix on its left (`A * v`) whereas a row
vector can be multiplied by a matrix on its right (such that `v.' * A = (A.' * v).'`). It
differs from a `1×n`-sized matrix by the facts that its transpose returns a vector and the
inner product `v1.' * v2` returns a scalar, but will otherwise behave similarly.
"""
struct RowVector{T,V<:AbstractVector} <: AbstractMatrix{T}
    vec::V
    function RowVector{T,V}(v::V) where V<:AbstractVector where T
        new(v)
    end
end

# Constructors that take a vector
@inline RowVector(vec::AbstractVector{T}) where {T} = RowVector{T,typeof(vec)}(vec)
@inline RowVector{T}(vec::AbstractVector{T}) where {T} = RowVector{T,typeof(vec)}(vec)

# Constructors that take a size and default to Array
@inline RowVector{T}(n::Int) where {T} = RowVector{T}(Vector{T}(n))
@inline RowVector{T}(n1::Int, n2::Int) where {T} = n1 == 1 ?
    RowVector{T}(Vector{T}(n2)) :
    error("RowVector expects 1×N size, got ($n1,$n2)")
@inline RowVector{T}(n::Tuple{Int}) where {T} = RowVector{T}(Vector{T}(n[1]))
@inline RowVector{T}(n::Tuple{Int,Int}) where {T} = n[1] == 1 ?
    RowVector{T}(Vector{T}(n[2])) :
    error("RowVector expects 1×N size, got $n")

# Conversion of underlying storage
convert(::Type{RowVector{T,V}}, rowvec::RowVector) where {T,V<:AbstractVector} =
    RowVector{T,V}(convert(V,rowvec.vec))

# similar tries to maintain the RowVector wrapper and the parent type
@inline similar(rowvec::RowVector) = RowVector(similar(parent(rowvec)))
@inline similar(rowvec::RowVector, ::Type{T}) where {T} = RowVector(similar(parent(rowvec), transpose_type(T)))

# Resizing similar currently loses its RowVector property.
@inline similar(rowvec::RowVector, ::Type{T}, dims::Dims{N}) where {T,N} = similar(parent(rowvec), T, dims)

parent(rowvec::RowVector) = rowvec.vec

# AbstractArray interface
@inline length(rowvec::RowVector) =  length(rowvec.vec)
@inline size(rowvec::RowVector) = (1, length(rowvec.vec))
@inline size(rowvec::RowVector, d) = ifelse(d==2, length(rowvec.vec), 1)
@inline indices(rowvec::RowVector) = (Base.OneTo(1), indices(rowvec.vec)[1])
@inline indices(rowvec::RowVector, d) = ifelse(d == 2, indices(rowvec.vec)[1], Base.OneTo(1))
IndexStyle(::RowVector) = IndexLinear()
IndexStyle(::Type{<:RowVector}) = IndexLinear()


@propagate_inbounds getindex(rowvec::RowVector, i) = rowvec.vec[i]
@propagate_inbounds setindex!(rowvec::RowVector, v, i) = setindex!(rowvec.vec, v, i)

# Cartesian indexing is distorted by getindex
# Furthermore, Cartesian indexes don't have to match shape, apparently!
@inline function getindex(rowvec::RowVector, i::CartesianIndex)
    @boundscheck if !(i.I[1] == 1 && i.I[2] ∈ indices(rowvec.vec)[1] && check_tail_indices(i.I...))
        throw(BoundsError(rowvec, i.I))
    end
    @inbounds return rowvec.vec[i.I[2]]
end
@inline function setindex!(rowvec::RowVector, v, i::CartesianIndex)
    @boundscheck if !(i.I[1] == 1 && i.I[2] ∈ indices(rowvec.vec)[1] && check_tail_indices(i.I...))
        throw(BoundsError(rowvec, i.I))
    end
    @inbounds rowvec.vec[i.I[2]] = v
end

@propagate_inbounds getindex(rowvec::RowVector, ::CartesianIndex{0}) = getindex(rowvec)
@propagate_inbounds getindex(rowvec::RowVector, i::CartesianIndex{1}) = getindex(rowvec, i.I[1])

@propagate_inbounds setindex!(rowvec::RowVector, v, ::CartesianIndex{0}) = setindex!(rowvec, v)
@propagate_inbounds setindex!(rowvec::RowVector, v, i::CartesianIndex{1}) = setindex!(rowvec, v, i.I[1])

# helper function for below
@inline to_vec(rowvec::RowVector) = rowvec.vec

# map: Preserve the RowVector by un-wrapping and re-wrapping, but note that `f`
# expects to operate within the transposed domain, so to_vec transposes the elements
@inline map(f, rowvecs::RowVector...) = RowVector(map(f, to_vecs(rowvecs...)...))

# broacast (other combinations default to higher-dimensional array)
@inline broadcast(f, rowvecs::Union{Number,RowVector}...) =
    RowVector(broadcast(f, to_vecs(rowvecs...)...))

# Horizontal concatenation #

@inline hcat(X::RowVector...) = RowVector(vcat(X...))
@inline hcat(X::Union{RowVector,Number}...) = RowVector(vcat(X...))

@inline typed_hcat(::Type{T}, X::RowVector...) where {T} =
    RowVector(typed_vcat(T, X...))
@inline typed_hcat(::Type{T}, X::Union{RowVector,Number}...) where {T} =
    RowVector(typed_vcat(T, X...))

# Multiplication #

# inner product -> dot product specializations
@inline *(rowvec::RowVector{T}, vec::AbstractVector{T}) where {T<:Real} = dotu(parent(rowvec), vec)

# Generic behavior
@inline function *(rowvec::RowVector, vec::AbstractVector)
    if length(rowvec) != length(vec)
        throw(DimensionMismatch("A has dimensions $(size(rowvec)) but B has dimensions $(size(vec))"))
    end
    sum(@inbounds(return rowvec[i]*vec[i]) for i = 1:length(vec))
end
@inline *(rowvec::RowVector, mat::AbstractMatrix) = RowVector(mat.' * rowvec.vec)
*(::RowVector, ::RowVector) = throw(DimensionMismatch("Cannot multiply two transposed vectors"))
@inline *(vec::AbstractVector, rowvec::RowVector) = vec .* rowvec



## Removed A_* overrides for now
