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
ArraySpace{SS<:Space,N}(sp::Array{SS,N}) =
    ArraySpace{SS,N,domaintype(first(sp)),mapreduce(rangetype,promote_type,sp)}(sp)
ArraySpace{N}(S::Space,n::NTuple{N,Int}) = ArraySpace(fill(S,n))
ArraySpace(S::Space,n::Integer) = ArraySpace(S,(n,))
ArraySpace(S::Space,n,m) = ArraySpace(fill(S,(n,m)))
ArraySpace(d::Domain,n...) = ArraySpace(Space(d),n...)

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

isambiguous(AS::ArraySpace) = isambiguous(AS.spaces[1])
## transforms


points(d::ArraySpace,n) = points(d.spaces[1],n)


transform{SS,V}(AS::ArraySpace{SS,1},vals::Vector{Vector{V}}) =
    transform(AS,transpose(hcat(vals...)))

#TODO: rework for different spaces
function transform{SS,T,V<:Number}(AS::ArraySpace{SS,1,T},M::Array{V,2})
    n=length(AS)

    @assert size(M,2)==n
    plan = plan_transform(AS.spaces[1],M[:,1])
    cfs=Vector{V}[plan*M[:,k]  for k=1:size(M,2)]

    interlace(cfs,AS)
end

# transform of array is same order as vectorizing and then transforming
transform{SS,n,V}(AS::ArraySpace{SS,n},vals::Vector{Array{V,n}}) =
    transform(vec(AS),map(vec,vals))
transform{SS,AV<:AbstractVector}(AS::ArraySpace{SS,1},vals::Vector{AV}) =
    transform(AS,map(Vector,vals))
transform{SS,n,V}(AS::ArraySpace{SS,1},vals::Vector{Vec{V,n}}) =
    transform(AS,map(Vector,vals))

Base.vec(AS::ArraySpace) = ArraySpace(vec(AS.spaces))
Base.vec{S,n,DD,RR}(f::Fun{ArraySpace{S,n,DD,RR}}) =
    [f[j] for j=1:length(f.space)]

Base.convert(::Type{Array},f::Fun{AS,T}) where {AS<:ArraySpace,T} =
    reshape(vec(f),size(space(f))...)

Base.map{AS<:ArraySpace}(f,A::Fun{AS}) = Base.collect_similar(A, Base.Generator(f,A))

Base.similar{SS,DD,RR}(a::Fun{ArraySpace{SS,1,DD,RR}}, S::Type) = Array{S,1}(size(a,1))
Base.similar{SS,DD,RR}(a::Fun{ArraySpace{SS,2,DD,RR}}, S::Type) = Array{S,2}(size(a,1), size(a,2))

Base.repmat(A::ArraySpace,n,m) = ArraySpace(repmat(A.spaces,n,m))

component(A::MatrixSpace,k::Integer,j::Integer) = A.spaces[k,j]

Base.getindex{DSS<:ArraySpace}(f::Fun{DSS},k::Integer) = component(f,k)


Base.getindex{S,DD,RR}(f::Fun{MatrixSpace{S,DD,RR}},k::Integer,j::Integer) =
    f[k+stride(f,2)*(j-1)]

Base.getindex{S,DD,RR}(f::Fun{MatrixSpace{S,DD,RR}},
                       k::Union{Integer,Range,Colon},
                       j::Union{Integer,Range,Colon}) =
    Fun(Array(f)[k,j])


function Base.vcat(vin::Fun...)
    #  remove tuple spaces
    v=Vector{Fun}(0)
    for f in vin
        if isa(space(f),VectorSpace)
            push!(v,vec(f)...)
        else
            push!(v,f)
        end
    end


    S = ArraySpace(space.(v))
    Fun(S,interlace(v,S))
end

Base.vcat(v::Union{Fun,Number}...) = vcat(map(Fun,v)...)

function Fun{F<:Fun}(v::Vector{F})
    S = ArraySpace(space.(v))
    Fun(S,interlace(v,S))
end

Space{S<:Space}(spl::Vector{S}) = ArraySpace(spl)


#TODO: rewrite
function Fun(v::Array{FF}) where {FF<:Fun}
    ff=Fun(vec(v))  # A vectorized version
    Fun(ArraySpace(map(space,v)),coefficients(ff))
end


function Fun{S,V,VV,DD,RR}(A::Array{Fun{VectorSpace{S,DD,RR},V,VV},2})
    @assert size(A,1)==1

    M=Matrix{Fun{S,V,VV}}(length(space(A[1])),size(A,2))
    for k=1:size(A,2)
        M[:,k]=vec(A[k])
    end
    Fun(M)
end

Fun(v::Array{NN}) where {NN<:Number} = Fun(v,ArraySpace(ConstantSpace(),size(v)...))

# Fun{SS,n}(v::Array{Any,n},sp::ArraySpace{SS,n}) = Fun(map((f,s)->Fun(f,s),v,sp))


# convert a vector to a Fun with ArraySpace

function Fun{TT,SS,n}(v::Array{TT,n},sp::ArraySpace{SS,n})
    if size(v) ≠ size(sp)
        throw(DimensionMismatch("Cannot convert $v to a Fun in space $sp"))
    end
    Fun(map(Fun,v,sp.spaces))
end
coefficients{TT,SS,n}(v::Array{TT,n},sp::ArraySpace{SS,n}) = coefficients(Fun(v,sp))


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

for OP in (:(Base.transpose),)
    @eval $OP{AS<:ArraySpace,T}(f::Fun{AS,T}) = Fun($OP(Array(f)))
end


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

coefficients(f::Vector,a::VectorSpace,b::VectorSpace) =
    interlace(map(coefficients,Fun(a,f),b),b)

coefficients{F<:Fun}(Q::Vector{F},rs::VectorSpace) =
    interlace(map(coefficients,Q,rs),rs)




Fun{FF<:Fun}(f::Vector{FF},d::VectorSpace) = Fun(d,coefficients(f,d))
Fun{FF<:Fun}(f::Matrix{FF},d::MatrixSpace) = Fun(d,coefficients(f,d))





## constructor

# change to ArraySpace
Fun{AS<:ArraySpace}(f::Fun{AS},d::ArraySpace) = space(f)==d ? f : Fun(d,coefficients(f,d))
Fun{AS<:ArraySpace}(f::Fun{AS},d::Space) = Fun(f,ArraySpace(d,size(space(f))))


# columns are coefficients
function Fun{T<:Number}(M::Array{T,2},sp::MatrixSpace)
    if size(M) ≠ size(sp)
        throw(DimensionMismatch())
    end
    Fun(map((f,s)->Fun(f,s),M,sp.spaces))
end

Fun(M::UniformScaling,sp::MatrixSpace) = Fun(M.λ*eye(size(sp)...),sp)


Fun{T<:Number}(M::Array{T,2},sp::Space) = Fun([Fun(M[:,k],sp) for k=1:size(M,2)])



Base.ones{T<:Number}(::Type{T},A::ArraySpace) = Fun(ones.(T,spaces(A)))
Base.ones(A::ArraySpace) = Fun(ones.(spaces(A)))

## calculus

for op in (:differentiate,:integrate,:(Base.cumsum),:(Base.real),:(Base.imag),:(Base.conj))
    @eval $op{V<:ArraySpace}(f::Fun{V}) = Fun(map($op,f))
end

function Base.det{A<:ArraySpace,V}(f::Fun{A,V})
    @assert size(space(f))==(2,2)
    m=Array(f)
    m[1,1]*m[2,2]-m[1,2]*m[2,1]
end

function Base.inv{A<:ArraySpace,T}(V::Fun{A,T})
    n,m = size(space(V))
    if n ≠ m
        throw(DimensionMismatch("space $(space(V)) is not square"))
    end

    # TODO: This assumes other columns have same spaces
    M=Multiplication(V,ArraySpace(space(V).spaces[:,1]))
    # convert I to the rangespace of M
    M\Fun(eye(m),repmat(rangespace(M),1,m))
end

## Algebra


const ArrayFun = Fun{S} where {S<:Space{D,R}} where {D,R<:AbstractArray}
const ScalarFun = Fun{S} where {S<:Space{D,R}} where {D,R<:Number}

for OP in (:*,:+,:-)
    @eval begin
        $OP(A::Array{<:Number},f::ArrayFun) = Fun($OP(A,Array(f)))
        $OP(f::ArrayFun,       A::Array{<:Number}) = Fun($OP(Array(f),A))
        $OP(A::Array{<:Fun},   f::ArrayFun) = Fun($OP(A,Array(f)))
        $OP(f::ArrayFun,       A::Array{<:Fun}) = Fun($OP(Array(f),A))
        $OP(A::UniformScaling, f::ArrayFun) = Fun($OP(A,Array(f)))
        $OP(f::ArrayFun,       A::UniformScaling) = Fun($OP(Array(f),A))
        $OP(A::Number,         f::ArrayFun) = Fun($OP(A,Array(f)))
        $OP(f::ArrayFun,       A::Number) = Fun($OP(Array(f),A))

        $OP(f::ScalarFun,      A::Array) = Fun(broadcast($OP,f,A))
        $OP(A::Array,          f::ScalarFun) = Fun(broadcast($OP,A,f))

        $OP(f::ScalarFun,      A::ArrayFun) = $OP(f,Array(A))
        $OP(A::ArrayFun,       f::ScalarFun) = $OP(Array(A),f)
    end
end

# use standard +, -
*(A::ArrayFun,f::ArrayFun) = Fun(Array(A)*Array(f))






# avoid ambiguity
# TODO: try removing this
for TYP in (:SpaceOperator,:TimesOperator,:QROperatorR,:QROperatorQ,:QROperator,:Operator)
    @eval \{S,DD,RR}(A::$TYP,b::Fun{MatrixSpace{S,DD,RR}};kwds...) =
        \(A,Array(b);kwds...)
end




## EuclideanSpace

const EuclideanSpace{RR} = VectorSpace{ConstantSpace{AnyDomain},AnyDomain,RR}
EuclideanSpace(n::Integer) = ArraySpace(ConstantSpace(),n)
