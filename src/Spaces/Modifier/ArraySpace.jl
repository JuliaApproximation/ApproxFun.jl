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

convert(::Type{Space},sp::AbstractArray{<:Space}) = ArraySpace(sp)
convert(::Type{Array},sp::ArraySpace) = sp.spaces
convert(::Type{Vector},sp::VectorSpace) = sp.spaces
convert(::Type{Matrix},sp::MatrixSpace) = sp.spaces


BlockInterlacer(sp::ArraySpace) = BlockInterlacer(blocklengths.(tuple(sp.spaces...)))
interlacer(sp::ArraySpace) = BlockInterlacer(sp)

for OP in (:(Base.length),:(Base.start),:(Base.endof),:(Base.size))
    @eval begin
        $OP(S::ArraySpace) = $OP(components(S))
        $OP(f::Fun{SS}) where {SS<:ArraySpace} = $OP(space(f))
    end
end

for OP in (:(getindex),:(Base.next),:(Base.done),:(Base.stride),:(Base.size))
    @eval $OP(S::ArraySpace,k) = $OP(components(S),k)
end


#support tuple set
for OP in (:(Base.done),:(Base.stride))
    @eval $OP(f::Fun{<:ArraySpace},k) = $OP(space(f),k)
end

getindex(f::ArraySpace,k...) = Space(component(f,k...))
Base.next(f::Fun{<:ArraySpace},k)=f[k],k+1


Base.reshape(AS::ArraySpace,k...) = ArraySpace(reshape(AS.spaces,k...))
dimension(AS::ArraySpace) = mapreduce(dimension,+,AS.spaces)

# TODO: union domain
domain(AS::ArraySpace) = domain(AS.spaces[1])
setdomain(A::ArraySpace,d::Domain) = ArraySpace(map(sp->setdomain(sp,d),A.spaces))



# support for Array of PiecewiseSpace


## transforms

#TODO: rework for different spaces
points(d::ArraySpace,n) = points(d.spaces[1],n)


transform(AS::ArraySpace{SS,1},vals::AbstractVector{Vector{V}}) where {SS,V} =
    transform(AS,transpose(hcat(vals...)))


function transform(AS::ArraySpace{SS,1,T},M::AbstractArray{V,2}) where {SS,T,V<:Number}
    n=length(AS)

    @assert size(M,2) == n
    plan = plan_transform(AS.spaces[1],M[:,1])
    cfs=Vector{V}[plan*M[:,k]  for k=1:size(M,2)]

    interlace(cfs,AS)
end

# transform of array is same order as vectorizing and then transforming
transform(AS::ArraySpace{SS,n},vals::AbstractVector{Array{V,n}}) where {SS,n,V} =
    transform(vec(AS),map(vec,vals))
transform(AS::VectorSpace{SS},vals::AbstractVector{AV}) where {SS,AV<:AbstractVector} =
    transform(AS,map(Vector,vals))
transform(AS::VectorSpace{SS},vals::AbstractVector{Vec{V,n}}) where {SS,n,V} =
    transform(AS,map(Vector,vals))

function itransform(AS::VectorSpace,cfs::AbstractVector)
    vf = vec(Fun(AS, cfs))
    n = maximum(ncoefficients, vf)
    vcat.(values.(pad!.(vf, n))...)
end


Base.vec(AS::ArraySpace) = ArraySpace(vec(AS.spaces))
Base.vec(f::Fun{ArraySpace{S,n,DD,RR}}) where {S,n,DD,RR} =
    [f[j] for j=1:length(f.space)]

Base.repmat(A::ArraySpace,n,m) = ArraySpace(repmat(A.spaces,n,m))

component(A::MatrixSpace,k::Integer,j::Integer) = A.spaces[k,j]

Base.getindex(f::Fun{DSS},k::Integer) where {DSS<:ArraySpace} = component(f,k)


Base.getindex(f::Fun{MatrixSpace{S,DD,RR}},k::Integer,j::Integer) where {S,DD,RR} =
    f[k+stride(f,2)*(j-1)]

Base.getindex(f::Fun{DSS},kj::CartesianIndex{1}) where {DSS<:ArraySpace} = f[kj[1]]
Base.getindex(f::Fun{DSS},kj::CartesianIndex{2}) where {DSS<:ArraySpace} = f[kj[1],kj[2]]


function Fun(A::AbstractArray{Fun{VectorSpace{S,DD,RR},V,VV},2}) where {S,V,VV,DD,RR}
    @assert size(A,1)==1

    M=Matrix{Fun{S,V,VV}}(length(space(A[1])),size(A,2))
    for k=1:size(A,2)
        M[:,k]=vec(A[k])
    end
    Fun(M)
end

# Fun{SS,n}(v::AbstractArray{Any,n},sp::ArraySpace{SS,n}) = Fun(map((f,s)->Fun(f,s),v,sp))


# convert a vector to a Fun with ArraySpace



function Fun(v::AbstractVector,sp::Space{D,R}) where {D,R<:AbstractVector}
    if size(v) ≠ size(sp)
        throw(DimensionMismatch("Cannot convert $v to a Fun in space $sp"))
    end
    Fun(map(Fun,v,components(sp)))
end

Fun(v::AbstractArray{TT,n},sp::Space{D,R}) where {D,R<:AbstractArray{SS,n}} where {TT,SS,n} =
    reshape(Fun(vec(v),vec(sp)),size(sp))


coefficients(v::AbstractArray{TT,n},sp::ArraySpace{SS,n}) where {TT,SS,n} = coefficients(Fun(v,sp))


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

function choosedomainspace(A::InterlaceOperator{T,1},sp::ArraySpace) where T
    # this ensures correct dispatch for unino
    sps = Vector{Space}(
        filter(x->!isambiguous(x),map(choosedomainspace,A.ops,sp.spaces)))
    if isempty(sps)
        UnsetSpace()
    else
        union(sps...)
    end
end


Base.reshape(f::Fun{AS},k...) where {AS<:ArraySpace} = Fun(reshape(space(f),k...),f.coefficients)

Base.diff(f::Fun{AS,T},n...) where {AS<:ArraySpace,T} = Fun(diff(Array(f),n...))

## conversion

function coefficients(f::AbstractVector,a::VectorSpace,b::VectorSpace)
    if size(a) ≠ size(b)
        throw(DimensionMismatch("dimensions must match"))
    end
    interlace(map(coefficients,Fun(a,f),b),b)
end


coefficients(Q::AbstractVector{F},rs::VectorSpace) where {F<:Fun} =
    interlace(map(coefficients,Q,rs),rs)




Fun(f::AbstractVector{FF},d::VectorSpace) where {FF<:Fun} = Fun(d,coefficients(f,d))
Fun(f::AbstractMatrix{FF},d::MatrixSpace) where {FF<:Fun} = Fun(d,coefficients(f,d))





## constructor



# columns are coefficients
function Fun(M::AbstractMatrix{<:Number},sp::MatrixSpace)
    if size(M) ≠ size(sp)
        throw(DimensionMismatch())
    end
    Fun(map((f,s)->Fun(f,s),M,sp.spaces))
end

Fun(M::UniformScaling,sp::MatrixSpace) = Fun(M.λ*eye(size(sp)...),sp)



Base.ones(::Type{T},A::ArraySpace) where {T<:Number} = Fun(ones.(T,spaces(A)))
Base.ones(A::ArraySpace) = Fun(ones.(spaces(A)))


## EuclideanSpace

const EuclideanSpace{RR} = VectorSpace{ConstantSpace{AnyDomain},AnyDomain,RR}
EuclideanSpace(n::Integer) = ArraySpace(ConstantSpace(Float64),n)




## support pieces

npieces(f::Fun{<:ArraySpace}) = npieces(f[1])
piece(f::Fun{<:ArraySpace}, k) = Fun(piece.(Array(f),k))
pieces(f::Fun{<:ArraySpace}) = [piece(f,k) for k=1:npieces(f)]
