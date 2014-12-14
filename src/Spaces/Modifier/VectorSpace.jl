## VectorSpace{T,S} encodes a space that is a Vector, with coefficients interlaced

export devec

immutable VectorDomainSpace{n,S,T} <: DomainSpace{T}
     space::S     
#      # for AnyDomain() usage
    VectorDomainSpace(sp::S)=new(sp)
    VectorDomainSpace(d::Domain)=new(S(d))
 end

VectorDomainSpace{T}(S::DomainSpace{T},n)=VectorDomainSpace{n,typeof(S),T}(S)
Base.length{n}(::VectorDomainSpace{n})=n

domain(S::VectorDomainSpace)=domain(S.space)
transform(S::VectorDomainSpace,vals::Vector)=transform!(S,hcat(vals...).')


function transform!{n}(S::VectorDomainSpace{n},M::Array)
    @assert size(M,2)==n
    for k=1:size(M,2)
        M[:,k]=transform(S.space,M[:,k])
    end
    vec(M.')
end

Base.vec{n,S<:DomainSpace,V,T}(f::Fun{VectorDomainSpace{n,S,V},T})=Fun{S,T}[Fun(f.coefficients[j:n:end],f.space.space) for j=1:n]

evaluate{V<:VectorDomainSpace,T}(f::Fun{V,T},x)=evaluate(vec(f),x)


# Base.ones{T<:Number,n}(::Type{T},S::VectorDomainSpace{n})=Fun(ones(T,n),S)









###########
# Piecewise Space
############

immutable PiecewiseSpace{S<:DomainSpace,T} <: DomainSpace{T}
    spaces::Vector{S} 
end
PiecewiseSpace(sp::Vector{Any})=PiecewiseSpace([sp...])
PiecewiseSpace{S,T}(::DomainSpace{T},spaces::Vector{S})=PiecewiseSpace{S,T}(spaces)
PiecewiseSpace(spaces)=PiecewiseSpace(first(spaces),spaces)
Space(d::UnionDomain)=PiecewiseSpace(map(Space,d.domains))
domain(S::PiecewiseSpace)=UnionDomain(map(domain,S.spaces))
Base.length(S::PiecewiseSpace)=S.spaces|>length
Base.getindex(d::PiecewiseSpace,k)=d.spaces[k]

Base.vec{S<:DomainSpace,V,T}(f::Fun{PiecewiseSpace{S,V},T},j::Integer)=Fun(f.coefficients[j:length(f.space):end],f.space.spaces[j])
Base.vec{S<:DomainSpace,V,T}(f::Fun{PiecewiseSpace{S,V},T})=Fun{S,T}[vec(f,j) for j=1:length(f.space)]




function spacescompatible{S,T}(A::PiecewiseSpace{S,T},B::PiecewiseSpace{S,T})
    if length(A) != length(B)
        false
    else
        ret=true
        for k=1:length(A)
            ret= ret && spacescompatible(A[k],B[k])
        end
        ret
    end
end




function transform{VV,ST,T}(S::PiecewiseSpace{VV,ST},vals::Vector{T})
    n=length(vals)
    K=length(S)
   k=div(n,K)
    PT=promote_type(ST,T)
    if k==0
        ret=Array(PT,n)
        for j=1:n
            ret[j]=transform(S[j],[vals[j]])[1]
        end
        
        ret
    else
        r=n-K*k
        M=Array(PT,k+1,K)
    
        for j=1:r
            M[:,j]=transform(S[j],vals[(j-1)*(k+1)+1:j*(k+1)])
        end
        for j=r+1:length(S)
            M[1:k,j]=transform(S[j],vals[r*(k+1)+(j-r-1)*k+1:r*(k+1)+(j-r)*k]) 
            M[k+1,j]=zero(PT)
        end    
        
    vec(M.')        
    end
end

itransform(S::PiecewiseSpace,cfs::Vector)=vcat([itransform(S.spaces[j],cfs[j:length(S):end]) for j=1:length(S)]...)


function evaluate{S<:PiecewiseSpace}(f::Fun{S},x::Number)
    d=domain(f)
    for k=1:length(d)
        if in(x,d[k])
            return vec(f,k)[x]
        end 
    end
end
evaluate{S<:PiecewiseSpace}(f::Fun{S},x::Vector)=[f[xk] for xk in x]

## space promotion

canonicalspace(sp::PiecewiseSpace)=PiecewiseSpace(map(canonicalspace,sp.spaces))

for op in (:maxspace,:minspace)
    @eval begin
        function $op(f::PiecewiseSpace,g::PiecewiseSpace)
            @assert length(f)==length(g)
            PiecewiseSpace([$op(f[k],g[k]) for k=1:length(f)])
        end
    end
end

for typ in (:PiecewiseSpace,:UnionDomain)
    @eval ==(a::($typ),b::($typ))=length(a)==length(b)&&all([a[k]==b[k] for k=1:length(a)])
end

## devec, asssume if domains the same we are vector




function devec{F<:Fun}(v::Vector{F})
    if spacescompatible(map(space,v))
        Fun(vec(coefficients(v).'),VectorDomainSpace(space(first(v)),length(v)))
    else
        Fun(vec(coefficients(v).'),PiecewiseSpace(map(space,v)))
    end
end

devec(v::Vector{Any})=devec([v...])

function devec{S<:FunctionSpace}(spl::Vector{S})
    if spacescompatible(spl)
        VectorDomainSpace(first(spl),length(spl))
    else
        PiecewiseSpace(spl)
    end
end

Base.vec(S::PiecewiseSpace)=S.spaces
Base.vec{n}(S::VectorDomainSpace{n})=fill(S.space,n)


for op in (:differentiate,:integrate)
    @eval $op{V<:Union(VectorDomainSpace,PiecewiseSpace)}(f::Fun{V})=devec(map($op,vec(f)))
end

Base.cumsum{V<:VectorDomainSpace}(f::Fun{V})=devec(map(cumsum,vec(f)))


function Base.cumsum{V<:PiecewiseSpace,T}(f::Fun{V,T})
    vf=vec(f)
    r=zero(T)
    for k=1:length(vf)
        vf[k]=cumsum(vf[k]) + r
        r=last(vf[k])
    end
    devec(vf)
end







## conversion



function spaceconversion{n}(f::Vector,a::VectorDomainSpace{n},b::VectorDomainSpace{n})
    A=a.space;B=b.space
    ret=copy(f)
    for k=1:n
        ret[k:n:end]=spaceconversion(ret[k:n:end],A,B)
    end
    ret
end

function spaceconversion{n}(f::Vector,a::VectorDomainSpace{n},b::PiecewiseSpace)
    A=a.space
    @assert n==length(b.spaces)
    ret=copy(f)
    for k=1:n
        ret[k:n:end]=spaceconversion(ret[k:n:end],A,b.spaces[k])
    end
    ret
end




## constructor

# columns are coefficients
Fun{T<:Number}(M::Array{T,2},sp::FunctionSpace)=devec([Fun(M[:,k],sp) for k=1:size(M,2)])


#There's no MatrixDomainSpace Yet
function Fun{T<:Number,n,S,Q}(M::Array{T,2},sp::VectorDomainSpace{n,S,Q})
    ret=Array(Fun{S,promote_type(Q,T)},n,size(M,2))
    for k=1:size(M,2)
        ret[:,k]=vec(Fun(M[:,k],sp))
    end
    ret
end




