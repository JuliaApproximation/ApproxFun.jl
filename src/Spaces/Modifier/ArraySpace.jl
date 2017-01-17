export devec,demat,mat


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
immutable ArraySpace{S,n,T,DD,dim} <: DirectSumSpace{NTuple{n,S},T,DD,dim}
     spaces::Array{S,n}
end

BlockInterlacer(sp::ArraySpace) = BlockInterlacer(blocklengths.(vec(sp.spaces)))
interlacer(sp::ArraySpace) = BlockInterlacer(sp)

typealias VectorSpace{S,T,DD,dim} ArraySpace{S,1,T,DD,dim}
typealias MatrixSpace{S,T,DD,dim} ArraySpace{S,2,T,DD,dim}

#TODO: Think through domain/domaindominsion
ArraySpace{SS<:Space,N}(sp::Array{SS,N}) =
    ArraySpace{SS,N,mapreduce(basistype,promote_type,sp),
               domaintype(first(sp)),domaindimension(first(sp))}(sp)
ArraySpace{N}(S::Space,n::NTuple{N,Int}) = ArraySpace(fill(S,n))
ArraySpace(S::Space,n::Integer) = ArraySpace(S,(n,))
ArraySpace(S::Space,n,m) = ArraySpace(fill(S,(n,m)))
ArraySpace(d::Domain,n...) = ArraySpace(Space(d),n...)



for FUNC in (:(Base.length),:(Base.size))
    @eval $FUNC(AS::ArraySpace) = $FUNC(AS.spaces)
end

for FUNC in (:(Base.size),:(Base.stride))
    @eval $FUNC(AS::ArraySpace,k) = $FUNC(AS.spaces,k)
end


Base.length{AS<:ArraySpace}(f::Fun{AS}) = length(space(f))
Base.stride{S,T,DD,dim}(AS::Fun{MatrixSpace{S,T,DD,dim}},k::Int) =
    k==1?k:size(AS,1)


Base.reshape(AS::ArraySpace,k...) = ArraySpace(reshape(AS.spaces,k...))
dimension(AS::ArraySpace) = mapreduce(dimension,*,AS.spaces)

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
Base.vec{S,n,T,DD,dim}(f::Fun{ArraySpace{S,n,T,DD,dim}}) =
    Fun{S,eltype(f)}[f[j] for j=1:length(f.space)]

mat{AS<:ArraySpace,T}(f::Fun{AS,T}) = reshape(vec(f),size(space(f))...)

# mat(f,1) vectorizes columnwise
function mat{S,V,T,DD,d}(f::Fun{MatrixSpace{S,V,DD,d},T},j::Integer)
    @assert j==1
    m=mat(f)
    r=Array(Fun{VectorSpace{S,V,DD,d},T},1,size(m,2))
    for k=1:size(m,2)
        r[1,k]=devec(m[:,k])
    end
    r
end


spaces(A::ArraySpace) = A.spaces
space(A::ArraySpace,k::Integer) = A.spaces[k]
space(A::MatrixSpace,k::Integer,j::Integer) = A.spaces[k,j]

Base.getindex{S,V,DD,d}(f::Fun{MatrixSpace{S,V,DD,d}},k::Integer,j::Integer) =
    f[k+stride(f,2)*(j-1)]

Base.getindex{S,V,DD,d}(f::Fun{MatrixSpace{S,V,DD,d}},k::Union{Integer,Range,Colon},j::Union{Integer,Range,Colon}) =
    Fun(mat(f)[k,j])

Base.getindex(S::ArraySpace,k::Integer) = S.spaces[k]
Base.getindex(S::ArraySpace,k::Integer,j::Integer) = S.spaces[k,j]

Base.start(S::ArraySpace) = start(S.spaces)
Base.next(S::ArraySpace,k) = next(S.spaces,k)
Base.done(S::ArraySpace,k) = done(S.spaces,k)
Base.endof(S::ArraySpace) = endof(S.spaces)


#support tuple set
for OP in (:(Base.start),:(Base.done),:(Base.endof))
    @eval $OP{SS<:ArraySpace}(f::Fun{SS},k...)=$OP(space(f),k...)
end

Base.next{SS<:ArraySpace}(f::Fun{SS},k)=f[k],k+1


function Base.vcat(vin::Fun...)
    #  remove tuple spaces
    v=Array(Fun,0)
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

function devec{F<:Fun}(v::Vector{F})
    S = ArraySpace(space.(v))
    Fun(S,interlace(v,S))
end

devec(v::Vector{Any}) = devec([v...])

devec{S<:Space}(spl::Vector{S}) = ArraySpace(spl)


#TODO: rewrite
function demat{FF<:Fun}(v::Array{FF})
    ff=devec(vec(v))  # A vectorized version
    Fun(ArraySpace(map(space,ff)),coefficients(ff))
end

demat(v::Vector{Any}) = devec(v)


function demat{S,T,V,DD,d}(A::Array{Fun{VectorSpace{S,T,DD,d},V},2})
    @assert size(A,1)==1

    M=Array(Fun{S,V},length(space(A[1])),size(A,2))
    for k=1:size(A,2)
        M[:,k]=vec(A[k])
    end
    demat(M)
end

Fun{F<:Fun}(V::AbstractVector{F}) = devec(V)
Fun{F<:Fun}(V::AbstractMatrix{F}) = demat(V)

Fun{SS,n}(v::Array{Any,n},sp::ArraySpace{SS,n}) = devec(map((f,s)->Fun(f,s),v,sp))


# convert a vector to a Fun with TupleSpace

function Fun{TT,SS,n}(v::Array{TT,n},sp::ArraySpace{SS,n})
    if size(v) ≠ size(sp)
        throw(DimensionMismatch("Cannot convert $v to a Fun in space $sp"))
    end
    demat(map(Fun,v,sp.spaces))
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
    @eval $OP{AS<:ArraySpace,T}(f::Fun{AS,T}) = demat($OP(mat(f)))
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

Base.diff{AS<:ArraySpace,T}(f::Fun{AS,T},n...) = demat(diff(mat(f),n...))

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
    demat(map((f,s)->Fun(f,s),M,sp.spaces))
end

Fun(M::UniformScaling,sp::MatrixSpace) = Fun(M.λ*eye(size(sp)...),sp)


Fun{T<:Number}(M::Array{T,2},sp::Space) = devec([Fun(M[:,k],sp) for k=1:size(M,2)])



Base.ones{T<:Number}(::Type{T},A::ArraySpace) = demat(ones.(T,spaces(A)))
Base.ones(A::ArraySpace) = demat(ones.(spaces(A)))

## calculus

for op in (:differentiate,:integrate,:(Base.cumsum),:(Base.real),:(Base.imag),:(Base.conj))
    @eval $op{V<:ArraySpace}(f::Fun{V}) = demat(map($op,f))
end

function Base.det{A<:ArraySpace,V}(f::Fun{A,V})
    @assert size(space(f))==(2,2)
    m=mat(f)
    m[1,1]*m[2,2]-m[1,2]*m[2,1]
end

function Base.inv{A<:ArraySpace,T}(V::Fun{A,T})
    if size(space(V),1) ≠ size(space(V),2)
        throw(DimensionMismatch("space $(space(V)) is not square"))
    end

    # TODO: This assumes other columns have same spaces
    M=Multiplication(V,ArraySpace(space(V).spaces[:,1]))
    # convert I to the rangespace of M
    M\Fun(eye(size(space(V),2)),ArraySpace(rangespace(M).space,size(space(V))))
end

## Algebra

for OP in (:*,:.*,:+,:-)
    @eval begin
        $OP{T<:Number,AS<:ArraySpace,V}(A::Array{T},f::Fun{AS,V}) = demat($OP(A,mat(f)))
        $OP{T<:Number,AS<:ArraySpace,V}(f::Fun{AS,V},A::Array{T}) = demat($OP(mat(f),A))
        $OP{T,S,AS<:ArraySpace,V}(A::Vector{Fun{S,T}},f::Fun{AS,V}) = demat($OP(A,mat(f)))
        $OP{T,S,AS<:ArraySpace,V}(f::Fun{AS,V},A::Vector{Fun{S,T}}) = demat($OP(mat(f),A))
        $OP{T,S,AS<:ArraySpace,V}(A::Array{Fun{S,T}},f::Fun{AS,V}) = demat($OP(A,mat(f)))
        $OP{T,S,AS<:ArraySpace,V}(f::Fun{AS,V},A::Array{Fun{S,T}}) = demat($OP(mat(f),A))
        $OP{AS<:ArraySpace,V}(A::UniformScaling,f::Fun{AS,V}) = demat($OP(A,mat(f)))
        $OP{AS<:ArraySpace,V}(f::Fun{AS,V},A::UniformScaling) = demat($OP(mat(f),A))
    end
end


for OP in (:*,:.*)
    @eval $OP{BS<:ArraySpace,T,AS<:ArraySpace,V}(A::Fun{BS,T},f::Fun{AS,V}) =
        demat($OP(mat(A),mat(f)))
end





## ConstantVectorSpace


typealias ConstantVectorSpace VectorSpace{ConstantSpace{AnyDomain},RealBasis,AnyDomain,1}


function Base.vec{V,TT,DD,d,T}(f::Fun{SumSpace{Tuple{ConstantVectorSpace,V},TT,DD,d},T},k)
    m=length(space(f)[1])
    if k≤m
        Fun(f.coefficients[k],ConstantSpace())
    else
        Fun(f.coefficients[m+1:end],space(f)[2])
    end
end



Base.vec{V,TT,DD,d,T}(f::Fun{SumSpace{Tuple{ConstantVectorSpace,V},TT,DD,d},T}) =
    Any[vec(f,k) for k=1:length(space(f)[1])+1]


# avoid ambiguity
for TYP in (:SpaceOperator,:TimesOperator,:QROperatorR,:QROperatorQ,:QROperator,:Operator)
    @eval \{S,T,DD,dim}(A::$TYP,b::Fun{MatrixSpace{S,T,DD,dim}};kwds...) =
        \(A,mat(b);kwds...)
end
