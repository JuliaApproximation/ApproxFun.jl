
export TensorSpace,⊗,ProductSpace

#  SV is a tuple of d spaces
abstract AbstractProductSpace{SV,T,d} <: Space{T,AnyDomain,d}


spacetype{SV}(::AbstractProductSpace{SV},k) = SV.parameters[k]


# Tensorizer
# This gives the map from coefficients to the
# tensor entry of a tensor product of d spaces
# findfirst is overriden to get efficient inverse
# blocklengths is a tuple of block lengths, e.g., Chebyshev()^2
# would be Tensorizer((1:∞,1:∞))
# ConstantSpace() ⊗ Chebyshev()
# would be Tensorizer((1:1,1:∞))
# and Chebyshev() ⊗ TupleSpace((Chebyshev(),Chebyshev()))
# would be Tensorizer((1:∞,2:2:∞))


immutable Tensorizer{DMS<:Tuple}
    blocks::DMS
end

Base.eltype(a::Tensorizer) = NTuple{length(a.blocks),Int}
Base.eltype{d,T}(::Tensorizer{NTuple{d,T}}) = NTuple{d,Int}
dimensions(a::Tensorizer) = map(sum,a.blocks)
Base.length(a::Tensorizer) = mapreduce(sum,*,a.blocks)

# (blockrow,blockcol), (subrow,subcol), (rowshift,colshift), (numblockrows,numblockcols), (itemssofar, length)
Base.start{DMS<:NTuple{2}}(a::Tensorizer{DMS}) = (1,1), (1,1), (0,0), (a.blocks[1][1],a.blocks[2][1]), (0,length(a))

function Base.next{DMS<:NTuple{2}}(a::Tensorizer{DMS},st)
    (K,J), (k,j), (rsh,csh), (n,m), (i,tot) = st
    ret = k+rsh,j+csh
    if k==n && j==m  # end of block
        if J == 1 || K == length(a.blocks[1])   # end of new block
            B = K+J # next block
            J = min(B, length(a.blocks[2]))  # don't go past new block
            K = B-J+1   # K+J-1 == B
        else
            K,J = K+1,J-1
        end
        k = j = 1
        if i+1 < tot # not done yet
            n,m = a.blocks[1][K], a.blocks[2][J]
            rsh,csh = sum(a.blocks[1][1:K-1]), sum(a.blocks[2][1:J-1])
        end
    elseif k==n
        k  = 1
        j += 1
    else
        k += 1
    end
    ret, ((K,J), (k,j), (rsh,csh), (n,m), (i+1,tot))
end


Base.done(a::Tensorizer,st) = st[end][1] ≥ st[end][2]

cache(a::Tensorizer) = CachedIterator(a)

function Base.findfirst(::Tensorizer{NTuple{2,Repeated{Bool}}},kj::Tuple{Int,Int})
    k,j=kj
    if k > 0 && j > 0
        n=k+j-2
        (n*(n+1))÷2+k
    else
        0
    end
end

# which block of the tensor
# equivalent to sum of indices -1

# block(it::Tensorizer,k) = sum(it[k])-length(it.blocks)+1
block{T}(ci::CachedIterator{T,Tensorizer{NTuple{2,Repeated{Bool}}}},k::Int) = sum(ci[k])-length(ci.iterator.blocks)+1

block(::Tensorizer{NTuple{2,Repeated{Bool}}},n::Int) =
    floor(Integer,sqrt(2n) + 1/2)
block(sp::Tensorizer,k::Int) = findfirst(x->x≥k,cumsum(blocklengths(sp)))
block(sp::CachedIterator,k::Int) = block(sp.iterator,k)

# [1,2,3] x 1:∞
function block(it::Tensorizer{Tuple{Vector{Bool},Repeated{Bool}}},n::Int)
    m=sum(it.blocks[1])
    if m == length(it.blocks[1])  # trivial blocks
        N=(m*(m+1))÷2
        if n < N
            return floor(Integer,sqrt(2n)+1/2)
        else
            return m+(n-N)÷m
        end
    else
        return findfirst(x->x≥n,cumsum(blocklengths(it)))
    end
end

# 1:∞ x 1:m
function block(it::Tensorizer{Tuple{Repeated{Bool},Vector{Bool}}},n::Int)
    m=length(it.blocks[2])  # assume all true
    N=(m*(m+1))÷2
    if n < N
        floor(Integer,sqrt(2n)+1/2)
    else
        m+(n-N)÷m
    end
end

blocklength(it,k) = blocklengths(it)[k]

blocklengths(::Tensorizer{NTuple{2,Repeated{Bool}}}) = 1:∞



blocklengths(it::Tensorizer) = tensorblocklengths(it.blocks...)
blocklengths(it::CachedIterator) = blocklengths(it.iterator)

function getindex(it::Tensorizer{NTuple{2,Repeated{Bool}}},n::Integer)
    m=block(it,n)
    p=findfirst(it,(1,m))
    j=1+n-p
    j,m-j+1
end





blockstart(it,K) = K==1?1:sum(blocklengths(it)[1:K-1])+1
blockstop(it,::Infinity{Bool}) = ∞
blockstop(it,K) = sum(blocklengths(it)[1:K])

blockstart(it,K::Block) = blockstart(it,K.K)
blockstop(it,K::Block) = blockstop(it,K.K)


blockrange(it,K) = blockstart(it,K):blockstop(it,K)




# convert from block, subblock to tensor
subblock2tensor(rt::Tensorizer{Tuple{Repeated{Bool},Repeated{Bool}}},K,k) =
    (k,K-k+1)

subblock2tensor{II}(rt::CachedIterator{II,Tensorizer{Tuple{Repeated{Bool},Repeated{Bool}}}},K,k) =
    (k,K-k+1)


subblock2tensor(rt::CachedIterator,K,k) = rt[blockstart(rt,K)+k-1]

# tensorblocklengths gives calculates the block sizes of each tensor product
#  Tensor product degrees are taken to be the sum of the degrees
#  a degree is which block you are in

tensorblocklengths(a) = a   # a single block is not modified
function tensorblocklengths(a::Repeated{Bool},b::Repeated{Bool})
    @assert a.x && b.x
    1:∞
end


function tensorblocklengths(a::Repeated,b::Repeated{Bool})
    @assert b.x
    a.x:a.x:∞
end


function tensorblocklengths(a::Repeated{Bool},b::Repeated)
    @assert a.x
    b.x:b.x:∞
end


function tensorblocklengths(a::Repeated,b::Repeated)
    m=a.x*b.x
    m:m:∞
end

function tensorblocklengths(a::Repeated,b)
    cs = a.x*cumsum(b)
    if isinf(length(b))
        cs
    elseif length(cs) == 1 && last(cs) == a.x
        a
    else
        flatten((cs,repeated(last(cs))))
    end
end



function tensorblocklengths(a::Repeated{Bool},b)
    @assert a.x
    cs = cumsum(b)
    if isinf(length(b))
        cs
    elseif length(cs) == 1 && last(cs) == a.x
        a
    else
        flatten((cs,repeated(last(cs))))
    end
end


tensorblocklengths(a,b::Repeated) =
    tensorblocklengths(b,a)


tensorblocklengths(a,b,c,d...) = tensorblocklengths(tensorblocklengths(a,b),c,d...)


# TensorSpace
# represents the tensor product of several subspaces

immutable TensorSpace{SV,T,d} <:AbstractProductSpace{SV,T,d}
    spaces::SV
end

tensorizer{SV,T,d}(sp::TensorSpace{SV,T,d}) = Tensorizer(map(blocklengths,sp.spaces))
blocklengths(S::TensorSpace) = tensorblocklengths(map(blocklengths,S.spaces)...)

TensorSpace(sp::Tuple) =
    TensorSpace{typeof(sp),mapreduce(basistype,promote_type,sp),
                mapreduce(domaindimension,+,sp)}(sp)


dimension(sp::TensorSpace) = mapreduce(dimension,*,sp.spaces)

for OP in (:spacescompatible,:(==))
    @eval $OP{SV,T,d}(A::TensorSpace{SV,T,d},B::TensorSpace{SV,T,d}) =
        all(Bool[$OP(A.spaces[k],B.spaces[k]) for k=1:length(A.spaces)])
end

canonicalspace(T::TensorSpace) = TensorSpace(map(canonicalspace,T.spaces))




coefficient_type(S::TensorSpace,T) =
    mapreduce(sp->coefficient_type(sp,T),promote_type,S.spaces)

TensorSpace(A...) = TensorSpace(tuple(A...))
TensorSpace(A::ProductDomain) = TensorSpace(tuple(map(Space,A.domains)...))
⊗(A::TensorSpace,B::TensorSpace) = TensorSpace(A.spaces...,B.spaces...)
⊗(A::TensorSpace,B::Space) = TensorSpace(A.spaces...,B)
⊗(A::Space,B::TensorSpace) = TensorSpace(A,B.spaces...)
⊗(A::Space,B::Space) = TensorSpace(A,B)

domain(f::TensorSpace) = mapreduce(domain,*,f.spaces)
Space(sp::ProductDomain) = TensorSpace(sp)

*(A::Space,B::Space) = A⊗B


# every column is in the same space for a TensorSpace
#TODO: remove
columnspace(S::TensorSpace,::) = S.spaces[1]

Base.length(d::TensorSpace) = length(d.spaces)
Base.getindex(d::TensorSpace,k::Integer) = d.spaces[k]


immutable ProductSpace{S<:Space,V<:Space,T} <: AbstractProductSpace{Tuple{S,V},T,2}
    spacesx::Vector{S}
    spacey::V
end

ProductSpace(spacesx::Vector,spacey)=ProductSpace{eltype(spacesx),
                                                  typeof(spacey),
                                                  promote_type(basistype(first(spacesx)),basistype(spacey))}(spacesx,spacey)

coefficient_type(S::ProductSpace,T) =
    promote_type(coefficient_type(S.spacesx[1],T),coefficient_type(S.spacesy,T))

⊗{S<:Space}(A::Vector{S},B::Space) = ProductSpace(A,B)
domain(f::ProductSpace) = domain(f.spacesx[1])*domain(f.spacesy)

Base.getindex(d::ProductSpace,k::Integer) = k==1?d.spacesx:d.spacey


space(d::AbstractProductSpace,k) = d[k]
isambiguous(A::TensorSpace) = isambiguous(A[1])||isambiguous(A[2])


Base.transpose(d::TensorSpace) = TensorSpace(d[2],d[1])





##Transforms

plan_column_transform(S,v) = plan_transform(columnspace(S,1),v)
plan_column_itransform(S,v) = plan_itransform(columnspace(S,1),v)

function itransform!(S::TensorSpace,M::Matrix)
    n=size(M,1)

    planc=plan_itransform(space(S,1),M[:,1])
    for k=1:size(M,2)
        M[:,k]=itransform(space(S,1),M[:,k],planc)
    end

    planr=plan_itransform(space(S,2),vec(M[1,:]))
    for k=1:n
        M[k,:]=itransform(space(S,2),vec(M[k,:]),planr)
    end
    M
end

function itransform!(S::AbstractProductSpace,M::Matrix)
    n=size(M,1)

    ## The order matters
    pln=plan_column_itransform(S,n)
    for k=1:size(M,2)
        M[:,k]=itransform(columnspace(S,k),M[:,k],pln)
    end

    for k=1:n
        M[k,:]=itransform(space(S,2),vec(M[k,:]))
    end
    M
end

function transform!(S::TensorSpace,M::Matrix)
    n=size(M,1)

    planc=plan_transform(space(S,1),M[:,1])
    for k=1:size(M,2)
        M[:,k]=transform(space(S,1),M[:,k],planc)
    end

    planr=plan_transform(space(S,2),vec(M[1,:]))
    for k=1:n
        M[k,:]=transform(space(S,2),vec(M[k,:]),planr)
    end
    M
end

function transform!{T}(S::AbstractProductSpace,M::Matrix{T})
    n=size(M,1)

    ## The order matters!!
    # For Disk Space, this is due to requiring decay
    # in function
    for k=1:n
        M[k,:]=transform(space(S,2),vec(M[k,:]))
    end

    pln=plan_column_transform(S,n)
    for k=1:size(M,2)
        # col may not be full length
        col=transform(columnspace(S,k),M[:,k],pln)
        M[1:length(col),k]=col
        for j=length(col)+1:n
            M[j,k]=zero(T) # fill rest with zeros
        end
    end


    M
end



## points

points(d::Union{BivariateDomain,BivariateSpace},n,m) =
    points(d,n,m,1),points(d,n,m,2)

function points(d::BivariateSpace,n,m,k)
    ptsx=points(columnspace(d,1),n)
    ptst=points(space(d,2),m)

    promote_type(eltype(ptsx),eltype(ptst))[fromcanonical(d,x,t)[k] for x in ptsx, t in ptst]
end




##  Fun routines

fromtensor(S::Space,M::Matrix) = fromtensor(tensorizer(S),M)
totensor(S::Space,M::Vector) = totensor(tensorizer(S),M)

# we only copy upper triangular of coefficients
function fromtensor(it::Tensorizer,M::Matrix)
    n,m=size(M)
    ret=zeros(eltype(M),blockstop(it,max(n,m)))
    k = 1
    for (K,J) in it
        if k > length(ret)
            break
        end
        if K ≤ n && J ≤ m
            ret[k] = M[K,J]
        end
        k += 1
    end
    ret
end


function totensor(it::Tensorizer,M::Vector)
    n=length(M)
    B=block(it,n)
    ds = dimensions(it)

    ret=zeros(eltype(M),sum(it.blocks[1][1:min(B,length(it.blocks[1]))]),
                        sum(it.blocks[2][1:min(B,length(it.blocks[2]))]))
    k=1
    for (K,J) in it
        if k > n
            break
        end
        ret[K,J] = M[k]
        k += 1
    end
    ret
end

for OP in (:block,:blockstart,:blockstop)
    @eval begin
        $OP(s::TensorSpace,::Infinity{Bool}) = ∞
        $OP(s::TensorSpace,M) = $OP(tensorizer(s),M)
    end
end

function points(sp::TensorSpace,n)
    pts=Array(Vec{2,Float64},0)
    if isfinite(dimension(sp[1])) && isfinite(dimension(sp[2]))
        N,M=dimension(sp[1]),dimension(sp[2])
    elseif isfinite(dimension(sp[1]))
        N=dimension(sp[1])
        M=n÷N
    elseif isfinite(dimension(sp[2]))
        M=dimension(sp[2])
        N=n÷M
    else
        N=M=round(Int,sqrt(n))
    end

    for y in points(sp[2],M),
        x in points(sp[1],N)
        push!(pts,Vec(x,y))
    end
    pts
end

function transform(sp::TensorSpace,vals,plan...)
    NM=length(vals)
    if isfinite(dimension(sp[1])) && isfinite(dimension(sp[2]))
        N,M=dimension(sp[1]),dimension(sp[2])
    elseif isfinite(dimension(sp[1]))
        N=dimension(sp[1])
        M=NM÷N
    elseif isfinite(dimension(sp[2]))
        M=dimension(sp[2])
        N=NM÷M
    else
        N=M=round(Int,sqrt(length(vals)))
    end

    V=reshape(copy(vals),N,M)

    fromtensor(sp,transform!(sp,V))
end

evaluate(f::AbstractVector,S::AbstractProductSpace,x) = ProductFun(totensor(S,f),S)(x...)
evaluate(f::AbstractVector,S::AbstractProductSpace,x,y) = ProductFun(totensor(S,f),S)(x,y)



coefficientmatrix{S<:AbstractProductSpace}(f::Fun{S}) = totensor(space(f),f.coefficients)



#TODO: Implement
# function ∂(d::TensorSpace{Interval{Float64}})
#     @assert length(d.spaces) ==2
#     PiecewiseSpace([d[1].a+im*d[2],d[1].b+im*d[2],d[1]+im*d[2].a,d[1]+im*d[2].b])
# end


union_rule(a::TensorSpace,b::TensorSpace) = TensorSpace(map(union,a.spaces,b.spaces))



## Convert from 1D to 2D


isconvertible(sp::UnivariateSpace,ts::TensorSpace) = length(ts.spaces) == 2 &&
    ((domain(ts)[1] == Point(0.0) && isconvertible(sp,ts[2])) ||
     (domain(ts)[2] == Point(0.0) && isconvertible(sp,ts[1])))


coefficients(f::Vector,sp::ConstantSpace,ts::TensorSpace) = f[1]*ones(ts).coefficients

function coefficients(f::Vector,sp::UnivariateSpace,ts::TensorSpace)
    @assert length(ts.spaces) == 2

    if domain(ts)[1] == Point(0.0)
        coefficients(f,sp,ts[2])
    elseif domain(ts)[2] == Point(0.0)
        coefficients(f,sp,ts[1])
    else
        error("Cannot convert coefficients from $sp to $ts")
    end
end


identity_fun(S::TensorSpace) = Fun(xyz->[xyz...],S)
