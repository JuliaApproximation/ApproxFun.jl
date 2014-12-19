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
    #TODO: Redesign
    if spacescompatible(spl)
        VectorDomainSpace(first(spl),length(spl))
    else
        PiecewiseSpace(spl)
    end
end

Base.vec{n}(S::VectorDomainSpace{n})=fill(S.space,n)


Base.cumsum{V<:VectorDomainSpace}(f::Fun{V})=devec(map(cumsum,vec(f)))






## conversion



function spaceconversion{n}(f::Vector,a::VectorDomainSpace{n},b::VectorDomainSpace{n})
    A=a.space;B=b.space
    ret=copy(f)
    for k=1:n
        ret[k:n:end]=spaceconversion(ret[k:n:end],A,B)
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




