# Creates a operator that permutes rows, in blocks of size
# length(perm)
struct PermutationOperator{T,DS,RS} <: Operator{T}
    perm::Vector{Int}
    domainspace::DS
    rangespace::RS
end
PermutationOperator{T}(prm, ds::DS, rs::RS) where {T, DS<:Space,RS<:Space} =
    PermutationOperator{T,DS,RS}(prm, ds, rs)
PermutationOperator(prm, ds, rs) = PermutationOperator{Int}(prm, ds, rs)
PermutationOperator(prm) = PermutationOperator(prm, ℓ⁰, ℓ⁰)

domainspace(P::PermutationOperator) = P.domainspace
rangespace(P::PermutationOperator) = P.rangespace


convert(::Type{Operator{T}},P::PermutationOperator) where {T} =
    PermutationOperator{T}(P.perm, P.domainspace, P.rangespace)

function bandwidths(P::PermutationOperator)
    dfs=P.perm-(1:length(P.perm))
    -minimum(dfs),maximum(dfs)
end

function getindex(P::PermutationOperator{T},k::Integer,j::Integer) where T
    n=length(P.perm)
    if (k-1)÷n == (j-1)÷n  # diagonal blocks
        k=mod(k-1,n)+1
        j=mod(j-1,n)+1
        convert(T,P.perm[k]==j)
    else
        zero(T)
    end
end

# the permutation that rearranges a to be b
multiplyperm(b,a) = Int[a[bk] for bk in b]
perm(a::Vector,b::Vector) = multiplyperm(invperm(sortperm(b)),sortperm(a))
perm(a::Tuple,b::Tuple) = perm(collect(a),collect(b))


struct NegateEven{T,DS,RS} <: Operator{T}
    domainspace::DS
    rangespace::RS
end

NegateEven{T}(ds::DS, rs::RS) where {DS<:Space,RS<:Space,T} =
    NegateEven{T,DS,RS}(ds, rs)
NegateEven(ds, rs) = NegateEven{Int}(ds, rs)
NegateEven() = NegateEven(ℓ⁰,ℓ⁰)

domainspace(P::NegateEven) = P.domainspace
rangespace(P::NegateEven) = P.rangespace

convert(::Type{Operator{T}},P::NegateEven) where {T} =
    NegateEven{T}(P.domainspace, P.rangespace)

bandwidths(P::NegateEven) = (0,0)


getindex(P::NegateEven{T},k::Integer,j::Integer) where {T} =
    k==j ? (iseven(k) ? -one(T) : one(T)) : zero(T)
