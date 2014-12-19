## VectorSpace{T,S} encodes a space that is a Vector, with coefficients interlaced

export devec

immutable VectorFunctionSpace{n,S,T,D<:Domain} <: FunctionSpace{T,D}
     space::S     
#      # for AnyDomain() usage
    VectorFunctionSpace(sp::S)=new(sp)
    VectorFunctionSpace(d::Domain)=new(S(d))
 end

VectorFunctionSpace{T,D}(S::FunctionSpace{T,D},n)=VectorFunctionSpace{n,typeof(S),T,D}(S)
Base.length{n}(::VectorFunctionSpace{n})=n

domain(S::VectorFunctionSpace)=domain(S.space)
transform(S::VectorFunctionSpace,vals::Vector)=transform!(S,hcat(vals...).')


function transform!{n}(S::VectorFunctionSpace{n},M::Array)
    @assert size(M,2)==n
    for k=1:size(M,2)
        M[:,k]=transform(S.space,M[:,k])
    end
    vec(M.')
end

Base.vec{n,S<:FunctionSpace,V,T,D<:Domain}(f::Fun{VectorFunctionSpace{n,S,V,D},T})=Fun{S,T}[Fun(f.coefficients[j:n:end],f.space.space) for j=1:n]

evaluate{V<:VectorFunctionSpace,T}(f::Fun{V,T},x)=evaluate(vec(f),x)


# Base.ones{T<:Number,n}(::Type{T},S::VectorFunctionSpace{n})=Fun(ones(T,n),S)



## devec, asssume if domains the same we are vector




function devec{F<:Fun}(v::Vector{F})
    if spacescompatible(map(space,v))
        Fun(vec(coefficients(v).'),VectorFunctionSpace(space(first(v)),length(v)))
    else
        Fun(vec(coefficients(v).'),PiecewiseSpace(map(space,v)))
    end
end

devec(v::Vector{Any})=devec([v...])

function devec{S<:FunctionSpace}(spl::Vector{S})
    #TODO: Redesign
    if spacescompatible(spl)
        VectorFunctionSpace(first(spl),length(spl))
    else
        PiecewiseSpace(spl)
    end
end

Base.vec{n}(S::VectorFunctionSpace{n})=fill(S.space,n)


Base.cumsum{V<:VectorFunctionSpace}(f::Fun{V})=devec(map(cumsum,vec(f)))






## conversion



function spaceconversion{n}(f::Vector,a::VectorFunctionSpace{n},b::VectorFunctionSpace{n})
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


#There's no MatrixFunctionSpace Yet
function Fun{T<:Number,n,S,Q}(M::Array{T,2},sp::VectorFunctionSpace{n,S,Q})
    ret=Array(Fun{S,promote_type(Q,T)},n,size(M,2))
    for k=1:size(M,2)
        ret[:,k]=vec(Fun(M[:,k],sp))
    end
    ret
end




