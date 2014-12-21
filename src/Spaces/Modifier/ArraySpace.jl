export devec,demat,mat



immutable ArraySpace{S,n,T,D<:Domain} <: FunctionSpace{T,D}
     space::S     
     dimensions::(Int...)
#      # for AnyDomain() usage
    ArraySpace(sp::S,dims)=new(sp,dims)
    ArraySpace(d::Domain,dims)=new(S(d),dims)
end


ArraySpace{T,D}(S::FunctionSpace{T,D},n::(Int...))=ArraySpace{typeof(S),length(n),T,D}(S,n)
ArraySpace{T,D}(S::FunctionSpace{T,D},n::Integer)=ArraySpace(S,(n,))
ArraySpace{T,D}(S::FunctionSpace{T,D},n,m)=ArraySpace{typeof(S),2,T,D}(S,(n,m))
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
    C=Array(promote_type(T,V),length(M))
    
    for k=1:size(M,2)
        C[k:n:end]=transform(AS.space,M[:,k])
    end
    C
end

# transform of array is same order as vectorizing and then transforming
transform{SS,n,V}(AS::ArraySpace{SS,n},vals::Vector{Array{V,n}})=transform(vec(AS),map(vec,vals))

Base.vec(AS::ArraySpace)=ArraySpace(AS.space,length(AS))
function Base.vec{S<:FunctionSpace,V,T,D<:Domain}(f::Fun{ArraySpace{S,1,V,D},T})
    n=length(space(f))
    Fun{S,T}[Fun(f.coefficients[j:n:end],space(f).space) for j=1:n]
end
Base.vec{AS<:ArraySpace,T}(f::Fun{AS,T})=vec(Fun(f.coefficients,vec(space(f))))

mat{AS<:ArraySpace,T}(f::Fun{AS,T})=reshape(vec(f),size(space(f))...)




function devec{F<:Fun}(v::Vector{F})
    @assert spacescompatible(map(space,v))
    
    Fun(vec(coefficients(v).'),ArraySpace(space(first(v)),length(v)))
end

devec(v::Vector{Any})=devec([v...])

function devec{S<:FunctionSpace}(spl::Vector{S})
    #TODO: Redesign
    @assert spacescompatible(spl)
    ArraySpace(first(spl),length(spl))
end


function demat{F<:Fun}(v::Array{F})
    ff=devec(vec(v))  # A vectorized version
    Fun(coefficients(ff),ArraySpace(space(ff).space,size(v)...))
end




## routines

spacescompatible(AS::ArraySpace,BS::ArraySpace)=size(AS)==size(BS) && spacescompatible(AS.space,BS.space)
canonicalspace(AS::ArraySpace)=ArraySpace(canonicalspace(AS.space),size(AS))
evaluate{AS<:ArraySpace,T}(f::Fun{AS,T},x)=evaluate(mat(f),x)
Base.transpose{AS<:ArraySpace,T}(f::Fun{AS,T})=demat(mat(f).')


## conversion

function spaceconversion{S,V,T}(f::Vector{T},a::ArraySpace{S,1},b::ArraySpace{V,1})
    n=length(a)
    @assert n==length(b)
    A=a.space;B=b.space
    ret=Array(T,length(f))
    for k=1:n
        ret[k:n:end]=spaceconversion(f[k:n:end],A,B)
    end
    ret
end






## constructor

# columns are coefficients
Fun{T<:Number}(M::Array{T,2},sp::FunctionSpace)=devec([Fun(M[:,k],sp) for k=1:size(M,2)])

# Automatically change to ArraySpace
Fun{T<:Number,S}(M::Array{T,2},sp::ArraySpace{S,1})=Fun(vec(M.'),ArraySpace(sp.space,length(sp),size(M,2))).'

## calculus

for op in (:differentiate,:integrate,:(Base.cumsum))
    @eval $op{V<:ArraySpace}(f::Fun{V})=demat(map($op,mat(f)))
end






## Algebra

*{T<:Number,AS<:ArraySpace,V}(A::Matrix{T},f::Fun{AS,V})=demat(A*mat(f))
*{BS<:ArraySpace,T,AS<:ArraySpace,V}(A::Fun{BS,T},f::Fun{AS,V})=demat(mat(A)*mat(f))
.*{BS<:ArraySpace,T,AS<:ArraySpace,V}(A::Fun{BS,T},f::Fun{AS,V})=demat(mat(A).*mat(f))


