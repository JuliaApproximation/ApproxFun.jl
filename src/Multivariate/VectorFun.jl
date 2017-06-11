const ArrayFun = Fun{S} where {S<:Space{D,R}} where {D,R<:AbstractArray}
const MatrixFun = Fun{S} where {S<:Space{D,R}} where {D,R<:AbstractMatrix}
const VectorFun = Fun{S} where {S<:Space{D,R}} where {D,R<:AbstractVector}
const ScalarFun = Fun{S} where {S<:Space{D,R}} where {D,R<:Number}


Base.convert(::Type{Array},f::ArrayFun) = reshape(vec(f),size(space(f))...)

Base.map(f,A::ArrayFun) = Base.collect_similar(A, Base.Generator(f,A))

Base.similar(a::VectorFun, S::Type) = Array{S,1}(size(a,1))
Base.similar(a::MatrixFun, S::Type) = Array{S,2}(size(a,1), size(a,2))


Base.getindex(f::MatrixFun,
               k::Union{Integer,Range,Colon},
               j::Union{Integer,Range,Colon}) =
    Fun(Array(f)[k,j])



function Base.vcat(vin::Fun...)
    #  remove tuple spaces
    v=Vector{Fun}(0)
    for f in vin
        if rangetype(space(f)) <: AbstractVector
            push!(v,vec(f)...)
        else
            push!(v,f)
        end
    end


    S = Space(space.(v))
    Fun(S,interlace(v,S))
end

Base.vcat(v::Union{Fun,Number}...) = vcat(map(Fun,v)...)

function Fun{F<:Fun}(v::AbstractVector{F})
    S = Space(space.(v))
    Fun(S,interlace(v,S))
end


#TODO: rewrite
function Fun(v::AbstractArray{<:Fun})
    ff=Fun(vec(v))  # A vectorized version
    Fun(Space(map(space,v)),coefficients(ff))
end

Fun(v::AbstractArray{NN}) where {NN<:Number} =
    Fun(v,Space(fill(ConstantSpace(NN),size(v))))
Fun(v::AbstractArray{Any}) = Fun(Fun.([v...]) :: AbstractArray{<:Fun})

Fun(f::ArrayFun,d::Space{D,R}) where {D,R<:AbstractArray} = space(f)==d ? f : Fun(d,coefficients(f,d))
Fun(f::ArrayFun,d::Space) = Fun(f,Space(fill(d,size(space(f)))))

Fun(M::AbstractMatrix{<:Number},sp::Space) = Fun([Fun(M[:,k],sp) for k=1:size(M,2)])

for OP in (:(Base.transpose),)
    @eval $OP(f::ArrayFun) = Fun($OP(Array(f)))
end

## calculus

for op in (:differentiate,:integrate,:(Base.cumsum),:(Base.real),:(Base.imag),:(Base.conj))
    @eval $op(f::ArrayFun) = Fun(map($op,f))
end

# TODO: use QR
function Base.det(f::MatrixFun)
    @assert size(space(f))==(2,2)
    m=Array(f)
    m[1,1]*m[2,2]-m[1,2]*m[2,1]
end

function Base.inv(V::MatrixFun)
    n,m = size(space(V))
    if n ≠ m
        throw(DimensionMismatch("space $(space(V)) is not square"))
    end

    # TODO: This assumes other columns have same spaces
    M=Multiplication(V,Space(space(V).spaces[:,1]))
    # convert I to the rangespace of M
    M\Fun(eye(m),repmat(rangespace(M),1,m))
end

## Algebra




for OP in (:*,:+,:-)
    @eval begin
        $OP(A::AbstractArray{<:Number}, f::ArrayFun) = Fun($OP(A,Array(f)))
        $OP(f::ArrayFun,                A::AbstractArray{<:Number}) = Fun($OP(Array(f),A))
        $OP(A::AbstractArray{<:Fun},    f::ArrayFun) = Fun($OP(A,Array(f)))
        $OP(f::ArrayFun,                A::AbstractArray{<:Fun}) = Fun($OP(Array(f),A))
        $OP(A::UniformScaling,          f::ArrayFun) = Fun($OP(A,Array(f)))
        $OP(f::ArrayFun,                A::UniformScaling) = Fun($OP(Array(f),A))
        $OP(A::Number,                  f::ArrayFun) = Fun($OP(A,Array(f)))
        $OP(f::ArrayFun,                A::Number) = Fun($OP(Array(f),A))

        $OP(f::ScalarFun,               A::AbstractArray) = Fun(broadcast($OP,f,A))
        $OP(A::AbstractArray,           f::ScalarFun) = Fun(broadcast($OP,A,f))

        $OP(f::ScalarFun,               A::ArrayFun) = $OP(f,Array(A))
        $OP(A::ArrayFun,                f::ScalarFun) = $OP(Array(A),f)
    end
end

# use standard +, -
*(A::ArrayFun,f::ArrayFun) = Fun(Array(A)*Array(f))






## Vector of fun routines

function coefficientmatrix{N,F}(::Type{N},f::AbstractVector{F},o...)
    if isempty(f)
        return Matrix{N}(0,0)
    end

    n=mapreduce(ncoefficients,max,f)
    m=length(f)
    R=zeros(N,n,m)
    for k=1:m
        R[1:ncoefficients(f[k]),k]=coefficients(f[k],o...)
    end
    R
end


scalarorfuntype{S,T<:Number}(::Fun{S,T}) = T
scalarorfuntype{T<:Number}(::T) = T
scalarorfuntype{T<:Number}(b::AbstractVector{T}) = T
scalarorfuntype(b::AbstractVector{Any}) = promote_type(map(scalarorfuntype,b)...)
scalarorfuntype{F<:Fun}(b::AbstractVector{F}) = promote_type(map(scalarorfuntype,b)...)


coefficientmatrix{F<:Fun}(Q::AbstractVector{F},o...)=coefficientmatrix(scalarorfuntype(Q),Q,o...)
coefficientmatrix(Q::AbstractVector{Any})=(@assert isempty(Q); zeros(0,0))


function values(f::AbstractVector{Fun{D,N,VN}}) where {D,N,VN}
    n=mapreduce(ncoefficients,max,f)
    m=length(f)
    R=zeros(N,n,m)
    for k=1:m
        R[:,k] = values(pad(f[k],n))
    end
    R
end

function values(p::AbstractMatrix{Fun{D,T,VT}}) where {D,T,VT}
    @assert size(p)[1] == 1

   values(vec(p))
end







## evaluation


#TODO: fix for complex
evaluate{T<:Fun}(A::AbstractArray{T},x::Number) =
    typeof(first(A)(x))[Akj(x) for Akj in A]



## Algebra
#
# ## scalar fun times vector
# #TODO: This is probably more efficient than the current implementation, but it
# # is too hacky to keep for now
#
#  for op in (:*,:(Base.Ac_mul_B),:(Base.At_mul_B))
#      @eval begin
#          function ($op)(A::Array{T,2}, p::AbstractVector{Fun{D,V,VT}}) where {T<:Number,V<:Number,D,VT}
#              cfs=$op(A,coefficientmatrix(p).')
#              ret = Vector{VFun{D,promote_type(T,V)}}(size(cfs,1))
#              for i = 1:size(cfs,1)
#                  ret[i] = chop!(Fun(first(p).space,cfs[i,:]),eps())
#              end
#              ret
#          end
#
#          function ($op)(p::AbstractVector{Fun{D,T,VT}},A::Array{T,2}) where {T<:Number,D,VT}
#              cfs=$op(A,coefficientmatrix(p).')
#              ret = Vector{VFun{D,T}}(size(cfs,1))
#              for i = 1:size(cfs,1)
#                  ret[i] = chop!(Fun(first(p).space,cfs[i,:]),eps())
#              end
#              ret
#          end
#      end
#  end
