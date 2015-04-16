export devec,demat,mat



immutable ArraySpace{S,n,T,dim} <: FunctionSpace{T,dim}
     space::S
     dimensions::(Int...)
#      # for AnyDomain() usage
    ArraySpace(sp::S,dims)=new(sp,dims)
    ArraySpace(d::Domain,dims)=new(S(d),dims)
end


ArraySpace{T,d}(S::FunctionSpace{T,d},n::(Int...))=ArraySpace{typeof(S),length(n),T,d}(S,n)
ArraySpace(S::FunctionSpace,n::Integer)=ArraySpace(S,(n,))
ArraySpace{T,d}(S::FunctionSpace{T,d},n,m)=ArraySpace{typeof(S),2,T,d}(S,(n,m))
ArraySpace(d::Domain,n...)=ArraySpace(Space(d),n...)
Base.length{SS}(AS::ArraySpace{SS,1})=AS.dimensions[1]
Base.length{SS}(AS::ArraySpace{SS,2})=*(AS.dimensions...)
Base.size(AS::ArraySpace)=AS.dimensions
Base.size(AS::ArraySpace,k)=AS.dimensions[k]

domain(AS::ArraySpace)=domain(AS.space)


## transforms


transform{SS,V}(AS::ArraySpace{SS,1},vals::Vector{Vector{V}})=transform(AS,hcat(vals...).')

function transform{SS,T,V<:Number}(AS::ArraySpace{SS,1,T},M::Array{V,2})
    n=length(AS)

    @assert size(M,2)==n
    cfs=Vector{V}[transform(AS.space,M[:,k])  for k=1:size(M,2)]

    C=zeros(coefficient_type(T,V),mapreduce(length,max,cfs)*size(M,2))

    for k=1:size(M,2)
        C[k:n:end]=cfs[k]
    end
    C
end

# transform of array is same order as vectorizing and then transforming
transform{SS,n,V}(AS::ArraySpace{SS,n},vals::Vector{Array{V,n}})=transform(vec(AS),map(vec,vals))

Base.vec(AS::ArraySpace)=ArraySpace(AS.space,length(AS))
function Base.vec{S<:FunctionSpace,V,T,d}(f::Fun{ArraySpace{S,1,V,d},T})
    n=length(space(f))
    Fun{S,T}[Fun(f.coefficients[j:n:end],space(f).space) for j=1:n]
end
Base.vec{AS<:ArraySpace,T}(f::Fun{AS,T})=vec(Fun(f.coefficients,vec(space(f))))

mat{AS<:ArraySpace,T}(f::Fun{AS,T})=reshape(vec(f),size(space(f))...)

# mat(f,1) vectorizes columnwise
function mat{S,V,T,d}(f::Fun{ArraySpace{S,2,V,d},T},j::Integer)
    @assert j==1
    m=mat(f)
    r=Array(Fun{ArraySpace{S,1,V,d},T},1,size(m,2))
    for k=1:size(m,2)
        r[1,k]=devec(m[:,k])
    end
    r
end




function devec{F<:Fun}(v::Vector{F})
    sps=map(space,v)
    if spacescompatible(sps)
        Fun(vec(coefficients(v).'),ArraySpace(first(sps),length(v)))
    else
        Fun(vec(coefficients(v).'),PiecewiseSpace(sps))
    end
end

devec(v::Vector{Any})=devec([v...])

function devec{S<:FunctionSpace}(spl::Vector{S})
    #TODO: Redesign
    @assert spacescompatible(spl)
    ArraySpace(first(spl),length(spl))
end


function demat{S,T}(v::Array{Fun{S,T}})
    ff=devec(vec(v))  # A vectorized version
    Fun(coefficients(ff),ArraySpace(space(ff).space,size(v)...))
end

demat(v::Vector{Any})=devec(v)


function demat{S,T,V}(A::Array{Fun{ArraySpace{S,1,T},V},2})
    @assert size(A,1)==1

    M=Array(Fun{S,V},length(space(A[1])),size(A,2))
    for k=1:size(A,2)
        M[:,k]=vec(A[k])
    end
    demat(M)
end





## routines

spacescompatible(AS::ArraySpace,BS::ArraySpace)=size(AS)==size(BS) && spacescompatible(AS.space,BS.space)
canonicalspace(AS::ArraySpace)=ArraySpace(canonicalspace(AS.space),size(AS))
evaluate{AS<:ArraySpace,T}(f::Fun{AS,T},x)=evaluate(mat(f),x)

for OP in (:(Base.transpose),)
    @eval $OP{AS<:ArraySpace,T}(f::Fun{AS,T})=demat($OP(mat(f)))
end

Base.diff{AS<:ArraySpace,T}(f::Fun{AS,T},n...)=demat(diff(mat(f),n...))

## conversion

function coefficients{S,V,T}(f::Vector{T},a::ArraySpace{S,1},b::ArraySpace{V,1})
    n=length(a)
    @assert n==length(b)
    A=a.space;B=b.space
    ret=Array(T,length(f))
    for k=1:n
        ret[k:n:end]=coefficients(f[k:n:end],A,B)
    end
    ret
end






## constructor

# columns are coefficients
Fun{T<:Number}(M::Array{T,2},sp::FunctionSpace)=devec([Fun(M[:,k],sp) for k=1:size(M,2)])

# Automatically change to ArraySpace
# A is interpreted as coefficients
function Fun{T<:Number,S}(A::Array{T,2},sp::ArraySpace{S,1})
    n=length(sp)
    m=size(A,2)

    cfs=Array(T,m*n*(div(size(A,1),n)+1))
    for k=1:size(A,1),j=1:m
        cfs[m*n*div(k-1,n)+mod(k-1,n)+(j-1)*n+1]=A[k,j]
    end
    # pad with zeros
    for k=size(A,1)+1:n*(div(size(A,1),n)+1),j=1:m
        cfs[m*n*div(k-1,n)+mod(k-1,n)+(j-1)*n+1]=zero(T)
    end
    Fun(cfs,ArraySpace(sp.space,n,size(A,2)))
end

Base.ones{T<:Number}(::Type{T},A::ArraySpace)=demat(fill(ones(T,A.space),A.dimensions...))
Base.ones(A::ArraySpace)=demat(fill(ones(A.space),A.dimensions...))

## calculus

for op in (:differentiate,:integrate,:(Base.cumsum))
    @eval $op{V<:ArraySpace}(f::Fun{V})=demat(map($op,mat(f)))
end


function Base.det{A<:ArraySpace,V}(f::Fun{A,V})
    @assert size(space(f))==(2,2)
    m=mat(f)
    m[1,1]*m[2,2]-m[1,2]*m[2,1]
end

function Base.inv{A<:ArraySpace,T}(V::Fun{A,T})
    @assert size(space(V),1)==size(space(V),2)
    M=Multiplication(V,ArraySpace(space(V).space,size(space(V),1)))
    M\eye(size(space(V),2))
end

## Algebra

for OP in (:*,:.*,:+,:-)
    @eval begin
        $OP{T<:Number,AS<:ArraySpace,V}(A::Array{T},f::Fun{AS,V})=demat($OP(A,mat(f)))
        $OP{T<:Number,AS<:ArraySpace,V}(f::Fun{AS,V},A::Array{T})=demat($OP(mat(f),A))
        $OP{T,S,AS<:ArraySpace,V}(A::Vector{Fun{S,T}},f::Fun{AS,V})=demat($OP(A,mat(f)))
        $OP{T,S,AS<:ArraySpace,V}(f::Fun{AS,V},A::Vector{Fun{S,T}})=demat($OP(mat(f),A))
        $OP{T,S,AS<:ArraySpace,V}(A::Array{Fun{S,T}},f::Fun{AS,V})=demat($OP(A,mat(f)))
        $OP{T,S,AS<:ArraySpace,V}(f::Fun{AS,V},A::Array{Fun{S,T}})=demat($OP(mat(f),A))
        $OP{AS<:ArraySpace,V}(A::UniformScaling,f::Fun{AS,V})=demat($OP(A,mat(f)))
        $OP{AS<:ArraySpace,V}(f::Fun{AS,V},A::UniformScaling)=demat($OP(mat(f),A))
    end
end


for OP in (:*,:.*)
    @eval $OP{BS<:ArraySpace,T,AS<:ArraySpace,V}(A::Fun{BS,T},f::Fun{AS,V})=demat($OP(mat(A),mat(f)))
end





## linsolve
# special implementation to solve column by column
function linsolve{S,T,Q}(A::BandedOperator,b::Fun{ArraySpace{S,2,T,1},Q};kwds...)
    rs=rangespace(A)
    if isa(rs,ArraySpace) && size(rs)==size(space(b))
        linsolve(A,[b];kwds...)
    else
        linsolve(A,mat(b,1);kwds...)
    end
end


