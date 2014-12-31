export depiece,pieces

###########
# Piecewise Space
############

immutable PiecewiseSpace{S<:FunctionSpace,T} <: FunctionSpace{T,UnionDomain}
    spaces::Vector{S} 
    PiecewiseSpace(::AnyDomain)=new(S[S(AnyDomain())])
    PiecewiseSpace(sp::Vector{S})=new(sp)
end
PiecewiseSpace(sp::Vector{Any})=PiecewiseSpace([sp...])
PiecewiseSpace{S,T}(::FunctionSpace{T},spaces::Vector{S})=PiecewiseSpace{S,T}(spaces)
PiecewiseSpace(spaces)=PiecewiseSpace(first(spaces),spaces)
Space(d::UnionDomain)=PiecewiseSpace(map(Space,d.domains))
domain(S::PiecewiseSpace)=UnionDomain(map(domain,S.spaces))
Base.length(S::PiecewiseSpace)=length(S.spaces)
Base.getindex(d::PiecewiseSpace,k)=d.spaces[k]

Base.vec{S<:FunctionSpace,V,T}(f::Fun{PiecewiseSpace{S,V},T},j::Integer)=Fun(f.coefficients[j:length(f.space):end],f.space.spaces[j])
Base.vec{S<:FunctionSpace,V,T}(f::Fun{PiecewiseSpace{S,V},T})=[vec(f,j) for j=1:length(f.space)]
pieces{S<:PiecewiseSpace,T}(f::Fun{S,T})=vec(f)
depiece{F<:Fun}(v::Vector{F})=Fun(vec(coefficients(v).'),PiecewiseSpace(map(space,v)))
depiece(v::Vector{Any})=depiece([v...])


Base.ones{T<:Number,SS,V}(::Type{T},S::PiecewiseSpace{SS,V})=depiece(Fun{SS,T}[ones(Sk) for Sk in S.spaces])





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
            if domain(f)==domain(g)
                PiecewiseSpace([$op(f[k],g[k]) for k=1:length(f)])
            else
                NoSpace()
            end
        end
    end
end

for typ in (:PiecewiseSpace,:UnionDomain)
    @eval ==(a::($typ),b::($typ))=length(a)==length(b)&&all([a[k]==b[k] for k=1:length(a)])
end




Base.vec(S::PiecewiseSpace)=S.spaces



## cumsum

for op in (:differentiate,:integrate)
    @eval $op{V<:PiecewiseSpace}(f::Fun{V})=depiece(map($op,pieces(f)))
end

for op in (:(Base.sign),)
    @eval $op{V<:PiecewiseSpace,T<:Real}(f::Fun{V,T})=depiece(map($op,pieces(f)))
end

for op in (:(.^),)
    @eval $op{V<:PiecewiseSpace}(f::Fun{V},k::Integer)=depiece(map(fk->$op(fk,k),pieces(f)))
    @eval $op{V<:PiecewiseSpace}(f::Fun{V},k)=depiece(map(fk->$op(fk,k),pieces(f)))
end

function Base.cumsum{V<:PiecewiseSpace,T}(f::Fun{V,T})
    vf=pieces(f)
    r=zero(T)
    for k=1:length(vf)
        vf[k]=cumsum(vf[k]) + r
        r=last(vf[k])
    end
    depiece(vf)
end




## Conversion from Vector to Piecewise
#  Right now if a Vector fun has different spaces in each component we represenent
#  it by a piecewise fun, so this allows conversion.


function spaceconversion{AS}(f::Vector,a::ArraySpace{AS,1},b::PiecewiseSpace)
    A=a.space
    n=length(a)
    @assert n==length(b.spaces)
    ret=copy(f)
    for k=1:n
        ret[k:n:end]=spaceconversion(ret[k:n:end],A,b.spaces[k])
    end
    ret
end
