


convert(::Type{Array},f::ArrayFun) = reshape(vec(f),size(space(f))...)

map(f,A::ArrayFun) = Base.collect_similar(A, Base.Generator(f,A))

similar(a::VectorFun, S::Type) = Array{S,1}(undef, size(a,1))
similar(a::MatrixFun, S::Type) = Array{S,2}(undef, size(a,1), size(a,2))


getindex(f::MatrixFun,
               k::Union{Integer,AbstractRange,Colon},
               j::Union{Integer,AbstractRange,Colon}) =
    Fun(Array(f)[k,j])


const FunTypes = Union{Fun,Number}
const ScalarFunTypes = Union{ScalarFun,Number}
function vcat(vin::FunTypes...)
    #  remove tuple spaces
    v=Vector{Fun}(undef,0)
    for f in vin
        if rangetype(space(f)) <: AbstractVector
            append!(v,vec(f))
        else
            push!(v,f)
        end
    end


    S = Space(space.(v))
    Fun(S,interlace(v,S))
end


function hcat(v::ScalarFunTypes...)
    ff = vcat(v...)  # A vectorized version
    transpose(ff)
end

hvcat(rows::Tuple{Vararg{Int}},v::FunTypes...) = Fun(hvnocat(rows,v...))


function hcat(v::VectorFun...)
    N = length(v[1])
    M = length(v)

    V = Array{Fun}(undef, N,M)
    for J=1:M
        V[:,J] = vec(v[J])
    end
    Fun(V)
end



function Fun(v::AbstractVector{F}) where F<:Fun
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
Fun(v::AbstractArray{Any}) = Fun(Fun.(collect(v)) :: AbstractArray{<:Fun})

Fun(f::ArrayFun,d::Space{D,R}) where {D,R<:AbstractArray} = space(f)==d ? f : Fun(d,coefficients(f,d))
Fun(f::ArrayFun,d::Space) = Fun(f,Space(fill(d,size(space(f)))))

Fun(M::AbstractMatrix{<:Number},sp::Space) = Fun([Fun(M[:,k],sp) for k=1:size(M,2)])

for OP in (:(transpose),)
    @eval begin
        $OP(f::ArrayFun) = Fun($OP(Array(f)))
        $OP(sp::Space{D,R}) where {D,R<:AbstractArray} = Space($OP(Array(sp)))
    end
end

## calculus

for op in (:differentiate,:integrate,:(cumsum),:(real),:(imag),:(conj))
    @eval $op(f::ArrayFun) = Fun(map($op,f))
end

# TODO: use QR
function det(f::MatrixFun)
    @assert size(space(f))==(2,2)
    m=Array(f)
    m[1,1]*m[2,2]-m[1,2]*m[2,1]
end

function inv(V::MatrixFun)
    n,m = size(space(V))
    if n ≠ m
        throw(DimensionMismatch("space $(space(V)) is not square"))
    end

    # TODO: This assumes other columns have same spaces
    M=Multiplication(V,Space(space(V).spaces[:,1]))
    # convert I to the rangespace of M
    M\Fun(Matrix(I,m,m), repeat(rangespace(M),1,m))
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

norm(A::VectorFun,p::Number) = norm(norm.(Array(A)),p)





## Vector of fun routines

function coefficientmatrix(::Type{N},f::AbstractVector{F},o...) where {N,F}
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


scalarorfuntype(::Fun{S,T}) where {S,T<:Number} = T
scalarorfuntype(::T) where {T<:Number} = T
scalarorfuntype(b::AbstractVector{T}) where {T<:Number} = T
scalarorfuntype(b::AbstractVector{Any}) = promote_type(map(scalarorfuntype,b)...)
scalarorfuntype(b::AbstractVector{F}) where {F<:Fun} = promote_type(map(scalarorfuntype,b)...)


coefficientmatrix(Q::AbstractVector{F},o...) where {F<:Fun}=coefficientmatrix(scalarorfuntype(Q),Q,o...)
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
evaluate(A::AbstractArray{T},x::Number) where {T<:Fun} =
    typeof(first(A)(x))[Akj(x) for Akj in A]
