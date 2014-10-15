## VectorSpace{T,S} encodes a space that is a Vector, with coefficients interlaced


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

Base.vec{S<:DomainSpace,V,T}(f::Fun{VectorDomainSpace{S,V},T})=Fun{S,T}[Fun(f.coefficients[j:length(f.space):end],f.space.space) for j=1:length(f.space)]

evaluate{V<:VectorDomainSpace,T}(f::Fun{V,T},x)=evaluate(vec(f),x)


# Base.ones{T<:Number,n}(::Type{T},S::VectorDomainSpace{n})=Fun(ones(T,n),S)
# Base.ones{O}(S::UltrasphericalSpace{O})=Fun(ones(1),S)    



## Separate domains

immutable UnionDomain{D<:Domain} <:Domain
    domains::Vector{D}
end

∪(d1::UnionDomain,d2::UnionDomain)=UnionDomain([d1.domains,d2.domains])
∪(d1::Domain,d2::UnionDomain)=UnionDomain([d1,d2.domains])
∪(d1::UnionDomain,d2::Domain)=UnionDomain([d1.domains,d2])
∪(d1::Domain,d2::Domain)=UnionDomain([d1,d2])
Base.length(d::UnionDomain)=d.domains|>length
Base.getindex(d::UnionDomain,k)=d.domains[k]
for op in (:(Base.first),:(Base.last))
    @eval $op(d::UnionDomain)=d.domains|>$op|>$op
end

function points(d::UnionDomain,n)
   k=div(n,length(d))
    r=n-length(d)*k

    [vcat([points(d.domains[j],k+1) for j=1:r]...),
        vcat([points(d.domains[j],k) for j=r+1:length(d)]...)]
end



immutable PiecewiseSpace{S<:DomainSpace,T} <: DomainSpace{T}
    spaces::Vector{S} 
end

PiecewiseSpace{S,T}(::DomainSpace{T},spaces::Vector{S})=PiecewiseSpace{S,T}(spaces)
PiecewiseSpace(spaces)=PiecewiseSpace(first(spaces),spaces)
Space(d::UnionDomain)=PiecewiseSpace(map(Space,d.domains))
domain(S::PiecewiseSpace)=UnionDomain(map(domain,S.spaces))
Base.length(S::PiecewiseSpace)=S.spaces|>length
Base.getindex(d::PiecewiseSpace,k)=d.spaces[k]

Base.vec{S<:DomainSpace,V,T}(f::Fun{PiecewiseSpace{S,V},T})=Fun{S,T}[Fun(f.coefficients[j:length(f.space):end],f.space.spaces[j]) for j=1:length(f.space)]



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



function transform{T}(S::PiecewiseSpace,vals::Vector{T})
    n=length(vals)
    K=length(S)
   k=div(n,K)
    r=n-K*k
    M=Array(Float64,k+1,K)
    
    for j=1:r
        M[:,j]=transform(S.spaces[j],vals[(j-1)*(k+1)+1:j*(k+1)])
    end
    for j=r+1:length(S)
        M[1:k,j]=transform(S.spaces[j],vals[r*(k+1)+(j-r-1)*k+1:r*(k+1)+(j-r)*k]) 
        M[k+1,j]=zero(T)
    end    
    vec(M.')
end

itransform(S::PiecewiseSpace,cfs::Vector)=vcat([itransform(S.spaces[j],cfs[j:length(S):end]) for j=1:length(S)]...)


function evaluate{S<:PiecewiseSpace}(f::Fun{S},x::Number)
    d=domain(f)
    vf=vec(f)
    for k=1:length(d)
        if in(x,d[k])
            return vf[k][x]
        end 
    end
end



## devec, asssume if domains the same we are vector




function devec{S,T}(v::Vector{Fun{S,T}})
    if spacescompatible(map(space,v))
        Fun(vec(coefficients(v).'),VectorDomainSpace(space(first(v)),length(v)))
    else
        Fun(vec(coefficients(v).'),PiecewiseSpace(map(space,v)))
    end
end




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




