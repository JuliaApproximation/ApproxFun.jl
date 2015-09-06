# Creates a operator that permutes rows, in blocks of size
# length(perm)
immutable PermutationOperator{T} <: BandedOperator{T}
    perm::Vector{Int}
end
PermutationOperator(prm)=PermutationOperator{Int}(prm)

Base.convert{BT<:Operator}(::Type{BT},P::PermutationOperator)=PermutationOperator{eltype(BT)}(P.perm)

function bandinds(P::PermutationOperator)
    dfs=P.perm-[1:length(P.perm);]
    minimum(dfs),maximum(dfs)
end

function addentries!(P::PermutationOperator,A,kr::Range)
    n=length(P.perm)
    for k=kr
        m=mod(k-1,n)+1
        A[k,P.perm[m]-m+k]+=1
    end
    A
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
