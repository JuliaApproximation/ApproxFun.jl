# Creates a operator that permutes rows, in blocks of size
# length(perm)
immutable PermutationOperator{T} <: Operator{T}
    perm::Vector{Int}
end
PermutationOperator(prm)=PermutationOperator{Int}(prm)

for OP in (:domainspace,:rangespace)
    @eval $OP(T::PermutationOperator) = ℓ⁰
end

Base.convert{T}(::Type{Operator{T}},P::PermutationOperator) =
    PermutationOperator{T}(P.perm)

function bandinds(P::PermutationOperator)
    dfs=P.perm-(1:length(P.perm))
    minimum(dfs),maximum(dfs)
end

function getindex{T}(P::PermutationOperator{T},k::Integer,j::Integer)
    n=length(P.perm)
    if (k-1)÷n == (j-1)÷n  # diagonal blocks
        k=mod(k-1,n)+1
        j=mod(j-1,n)+1
        T(P.perm[k]==j)
    else
        zero(T)
    end
end


multiplyperm(b,a)=Int[a[bk] for bk in b]
perm(a::Vector,b::Vector)=multiplyperm(invperm(sortperm(b)),sortperm(a))
perm(a::Tuple,b::Tuple)=perm([a...],[b...])

# the operator that rearranges a to be b
function PermutationOperator{T}(::Type{T},a::Vector,b::Vector)
    @assert sort(a)==sort(b)
    PermutationOperator{T}(perm(a,b))
end
PermutationOperator{T}(::Type{T},a::Tuple,b::Tuple)=PermutationOperator(T,[a...],[b...])


immutable NegateEven{T} <: Operator{T} end

NegateEven() = NegateEven{Float64}()

for OP in (:domainspace,:rangespace)
    @eval $OP(T::NegateEven) = ℓ⁰
end

Base.convert{T}(::Type{Operator{T}},P::NegateEven) =
    NegateEven{T}()

bandinds(P::NegateEven) = (0,0)


getindex{T}(P::NegateEven{T},k::Integer,j::Integer) = k==j ? (iseven(k)?-one(T):one(T)) : zero(T)
