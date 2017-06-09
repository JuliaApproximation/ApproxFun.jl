doc"""
    ArraySpace(s::Space,dims...)

is used to represent array-valued expansions in a space `s`.  The
coefficients are of each entry are interlaced.

For example,
```julia
f = Fun(x->[exp(x),sin(x)],-1..1)
space(f) == ArraySpace(Chebyshev(),2)
```
"""
# TODO: support general vector types
struct ArraySpace{S,n,DD,RR} <: DirectSumSpace{NTuple{n,S},DD,Array{RR,n}}
     spaces::Array{S,n}
end

const VectorSpace{S,DD,RR} = ArraySpace{S,1,DD,RR}
const MatrixSpace{S,DD,RR} = ArraySpace{S,2,DD,RR}

#TODO: Think through domain/domaindominsion
ArraySpace(sp::AbstractArray{SS,N}) where {SS<:Space,N} =
    ArraySpace{SS,N,domaintype(first(sp)),mapreduce(rangetype,promote_type,sp)}(sp)
ArraySpace(S::Space,n::NTuple{N,Int}) where {N} = ArraySpace(fill(S,n))
ArraySpace(S::Space,n::Integer) = ArraySpace(S,(n,))
ArraySpace(S::Space,n,m) = ArraySpace(fill(S,(n,m)))
ArraySpace(d::Domain,n...) = ArraySpace(Space(d),n...)

Base.convert(::Type{Space},sp::AbstractArray{<:Space}) = ArraySpace(sp)

BlockInterlacer(sp::ArraySpace) = BlockInterlacer(blocklengths.(tuple(sp.spaces...)))
interlacer(sp::ArraySpace) = BlockInterlacer(sp)

for OP in (:(Base.length),:(Base.start),:(Base.endof),:(Base.size))
    @eval begin
        $OP(S::ArraySpace) = $OP(components(S))
        $OP{SS<:ArraySpace}(f::Fun{SS}) = $OP(space(f))
    end
end

for OP in (:(Base.getindex),:(Base.next),:(Base.done),:(Base.stride),:(Base.size))
    @eval $OP(S::ArraySpace,k) = $OP(components(S),k)
end


#support tuple set
for OP in (:(Base.done),:(Base.stride))
    @eval $OP{SS<:ArraySpace}(f::Fun{SS},k) = $OP(space(f),k)
end

getindex(f,k...) = component(f,k...)
Base.next{SS<:ArraySpace}(f::Fun{SS},k)=f[k],k+1


Base.reshape(AS::ArraySpace,k...) = ArraySpace(reshape(AS.spaces,k...))
dimension(AS::ArraySpace) = mapreduce(dimension,+,AS.spaces)

# TODO: union domain
domain(AS::ArraySpace) = domain(AS.spaces[1])
setdomain(A::ArraySpace,d::Domain) = ArraySpace(map(sp->setdomain(sp,d),A.spaces))


isambiguous(AS::ArraySpace) = isambiguous(AS.spaces[1])
## transforms


points(d::ArraySpace,n) = points(d.spaces[1],n)


transform{SS,V}(AS::ArraySpace{SS,1},vals::AbstractVector{Vector{V}}) =
    transform(AS,transpose(hcat(vals...)))

#TODO: rework for different spaces
function transform{SS,T,V<:Number}(AS::ArraySpace{SS,1,T},M::AbstractArray{V,2})
    n=length(AS)

    @assert size(M,2)==n
    plan = plan_transform(AS.spaces[1],M[:,1])
    cfs=Vector{V}[plan*M[:,k]  for k=1:size(M,2)]

    interlace(cfs,AS)
end

# transform of array is same order as vectorizing and then transforming
transform{SS,n,V}(AS::ArraySpace{SS,n},vals::AbstractVector{Array{V,n}}) =
    transform(vec(AS),map(vec,vals))
transform{SS,AV<:AbstractVector}(AS::ArraySpace{SS,1},vals::AbstractVector{AV}) =
    transform(AS,map(Vector,vals))
transform{SS,n,V}(AS::ArraySpace{SS,1},vals::AbstractVector{Vec{V,n}}) =
    transform(AS,map(Vector,vals))

Base.vec(AS::ArraySpace) = ArraySpace(vec(AS.spaces))
Base.vec{S,n,DD,RR}(f::Fun{ArraySpace{S,n,DD,RR}}) =
    [f[j] for j=1:length(f.space)]

Base.repmat(A::ArraySpace,n,m) = ArraySpace(repmat(A.spaces,n,m))

component(A::MatrixSpace,k::Integer,j::Integer) = A.spaces[k,j]

Base.getindex(f::Fun{DSS},k::Integer) where {DSS<:ArraySpace} = component(f,k)


Base.getindex{S,DD,RR}(f::Fun{MatrixSpace{S,DD,RR}},k::Integer,j::Integer) =
    f[k+stride(f,2)*(j-1)]

Base.getindex(f::Fun{DSS},kj::CartesianIndex{1}) where {DSS<:ArraySpace} = f[kj[1]]
Base.getindex(f::Fun{DSS},kj::CartesianIndex{2}) where {DSS<:ArraySpace} = f[kj[1],kj[2]]


function Fun{S,V,VV,DD,RR}(A::AbstractArray{Fun{VectorSpace{S,DD,RR},V,VV},2})
    @assert size(A,1)==1

    M=Matrix{Fun{S,V,VV}}(length(space(A[1])),size(A,2))
    for k=1:size(A,2)
        M[:,k]=vec(A[k])
    end
    Fun(M)
end

# Fun{SS,n}(v::AbstractArray{Any,n},sp::ArraySpace{SS,n}) = Fun(map((f,s)->Fun(f,s),v,sp))


# convert a vector to a Fun with ArraySpace

function Fun{TT,SS,n}(v::AbstractArray{TT,n},sp::ArraySpace{SS,n})
    if size(v) ≠ size(sp)
        throw(DimensionMismatch("Cannot convert $v to a Fun in space $sp"))
    end
    Fun(map(Fun,v,sp.spaces))
end
coefficients{TT,SS,n}(v::AbstractArray{TT,n},sp::ArraySpace{SS,n}) = coefficients(Fun(v,sp))


for (OPrule,OP) in ((:conversion_rule,:conversion_type),(:maxspace_rule,:maxspace),
                        (:union_rule,:union))
    # ArraySpace doesn't allow reordering
    @eval function $OPrule(S1::ArraySpace,S2::ArraySpace)
        sps = map($OP,S1.spaces,S2.spaces)
        for s in sps
            if isa(s,NoSpace)
                return NoSpace()
            end
        end
        ArraySpace(sps)
    end
end

## routines

spacescompatible(AS::ArraySpace,BS::ArraySpace) =
    size(AS) == size(BS) && all(spacescompatible.(AS.spaces,BS.spaces))
canonicalspace(AS::ArraySpace) = ArraySpace(canonicalspace.(AS.spaces))
evaluate(f::AbstractVector,S::ArraySpace,x) = map(g->g(x),Fun(S,f))


## choosedomainspace

function choosedomainspace{T}(A::InterlaceOperator{T,1},sp::ArraySpace)
    # this ensures correct dispatch for unino
    sps = Vector{Space}(
        filter(x->!isambiguous(x),map(choosedomainspace,A.ops,sp.spaces)))
    if isempty(sps)
        UnsetSpace()
    else
        union(sps...)
    end
end


Base.reshape{AS<:ArraySpace}(f::Fun{AS},k...) = Fun(reshape(space(f),k...),f.coefficients)

Base.diff{AS<:ArraySpace,T}(f::Fun{AS,T},n...) = Fun(diff(Array(f),n...))

## conversion

coefficients(f::AbstractVector,a::VectorSpace,b::VectorSpace) =
    interlace(map(coefficients,Fun(a,f),b),b)

coefficients{F<:Fun}(Q::AbstractVector{F},rs::VectorSpace) =
    interlace(map(coefficients,Q,rs),rs)




Fun{FF<:Fun}(f::AbstractVector{FF},d::VectorSpace) = Fun(d,coefficients(f,d))
Fun{FF<:Fun}(f::AbstractMatrix{FF},d::MatrixSpace) = Fun(d,coefficients(f,d))





## constructor



# columns are coefficients
function Fun(M::AbstractMatrix{<:Number},sp::MatrixSpace)
    if size(M) ≠ size(sp)
        throw(DimensionMismatch())
    end
    Fun(map((f,s)->Fun(f,s),M,sp.spaces))
end

Fun(M::UniformScaling,sp::MatrixSpace) = Fun(M.λ*eye(size(sp)...),sp)



Base.ones{T<:Number}(::Type{T},A::ArraySpace) = Fun(ones.(T,spaces(A)))
Base.ones(A::ArraySpace) = Fun(ones.(spaces(A)))


## EuclideanSpace

const EuclideanSpace{RR} = VectorSpace{ConstantSpace{AnyDomain},AnyDomain,RR}
EuclideanSpace(n::Integer) = ArraySpace(ConstantSpace(),n)
