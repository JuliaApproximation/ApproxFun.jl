export PiecewiseSpace,depiece,pieces

###########
# Piecewise Space
############

immutable PiecewiseSpace{S,T,d} <: Space{T,AnyDomain,d}
    spaces::Vector{S}
    PiecewiseSpace(::AnyDomain)=new(S[S(AnyDomain())])
    PiecewiseSpace(sp::Vector{S})=new(sp)
end
PiecewiseSpace(sp::Vector{Any})=PiecewiseSpace([sp...])
PiecewiseSpace(S::Space,spaces::Vector)=PiecewiseSpace{eltype(spaces),basistype(S),ndims(S)}(spaces)
PiecewiseSpace(spaces)=PiecewiseSpace(first(spaces),spaces)
PiecewiseSpace(a::Space,b::Space)=PiecewiseSpace([a,b])
Space(d::UnionDomain)=PiecewiseSpace(map(Space,d.domains))
domain(S::PiecewiseSpace)=UnionDomain(map(domain,S.spaces))
Base.length(S::PiecewiseSpace)=length(S.spaces)
Base.getindex(d::PiecewiseSpace,k)=d.spaces[k]

Base.vec{S<:Space,V,T,d}(f::Fun{PiecewiseSpace{S,V,d},T},j::Integer)=Fun(f.coefficients[j:length(f.space):end],f.space.spaces[j])
Base.vec{S<:Space,V,T,d}(f::Fun{PiecewiseSpace{S,V,d},T})=[vec(f,j) for j=1:length(f.space)]
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



function points(d::PiecewiseSpace,n)
   k=div(n,length(d))
    r=n-length(d)*k

    [vcat([points(d.spaces[j],k+1) for j=1:r]...);
        vcat([points(d.spaces[j],k) for j=r+1:length(d)]...)]
end



function transform(S::PiecewiseSpace,vals::Vector,plan...)
    n=length(vals)
    K=length(S)
    k=div(n,K)
    PT=coefficient_type(S,eltype(vals))
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

itransform(S::PiecewiseSpace,cfs::Vector,plan...)=vcat([itransform(S.spaces[j],cfs[j:length(S):end]) for j=1:length(S)]...)


function evaluate{S<:PiecewiseSpace}(f::Fun{S},x::Number)
    d=domain(f)
    for k=1:length(d)
        if in(x,d[k])
            return vec(f,k)(x)
        end
    end
end
evaluate{S<:PiecewiseSpace}(f::Fun{S},x::Vector)=[f(xk) for xk in x]

## space promotion

canonicalspace(sp::PiecewiseSpace)=PiecewiseSpace(map(canonicalspace,sp.spaces))

for op in (:maxspace,:conversion_type)
    @eval begin
        function $op(f::PiecewiseSpace,g::PiecewiseSpace)
            if domain(f)==domain(g)
                # hack to make sure type is correct
                PiecewiseSpace([[$op(f[k],g[k]) for k=1:length(f)]...])
            else
                NoSpace()
            end
        end
    end
end

for typ in (:PiecewiseSpace,:UnionDomain)
    @eval ==(a::($typ),b::($typ))=length(a)==length(b)&&all(k->a[k]==b[k],1:length(a))
end




Base.vec(S::PiecewiseSpace)=S.spaces



## cumsum

for op in (:differentiate,:integrate)
    @eval $op{V<:PiecewiseSpace}(f::Fun{V})=depiece(map($op,pieces(f)))
end



Base.dot{S<:PiecewiseSpace,V<:PiecewiseSpace}(f::Fun{S},g::Fun{V}) = sum(map(dot,pieces(f),pieces(g)))

function Base.cumsum{V<:PiecewiseSpace,T}(f::Fun{V,T})
    vf=pieces(f)
    r=zero(T)
    for k=1:length(vf)
        vf[k]=cumsum(vf[k]) + r
        r=last(vf[k])
    end
    depiece(vf)
end


Base.sum{V<:PiecewiseSpace,T}(f::Fun{V,T})=mapreduce(sum,+,pieces(f))



## Conversion from Vector to Piecewise
#  Right now if a Vector fun has different spaces in each component we represenent
#  it by a piecewise fun, so this allows conversion.


function coefficients(f::Vector,a::VectorSpace,b::PiecewiseSpace)
    A=a.space
    n=length(a)
    @assert n==length(b.spaces)
    ret=copy(f)
    for k=1:n
        ret[k:n:end]=coefficients(ret[k:n:end],A,b.spaces[k])
    end
    ret
end


union_rule(P::PiecewiseSpace,C::ConstantSpace)=PiecewiseSpace(map(sp->union(sp,C),P.spaces))


## Definite Integral

# This makes sure that the defaults from a given Domain are respected for the UnionDomain.

DefiniteIntegral(d::UnionDomain) = DefiniteIntegral(PiecewiseSpace(map(domainspace,map(DefiniteIntegral,d.domains))))
DefiniteLineIntegral(d::UnionDomain) = DefiniteLineIntegral(PiecewiseSpace(map(domainspace,map(DefiniteLineIntegral,d.domains))))

####### This is a hack to get the Faraday Cage working.
function getindex{PWS<:PiecewiseSpace,T}(Σ::DefiniteLineIntegral{PWS,T},kr::Range)
    d = domain(Σ)
    n = length(d)
    promote_type(T,eltype(d))[k ≤ n? one(T) : zero(T) for k=kr]
end
datalength{PWS<:PiecewiseSpace,T}(Σ::DefiniteLineIntegral{PWS,T})=length(domain(Σ))
####### This is a hack to get the Faraday Cage working.

## TensorSpace of two PiecewiseSpaces

Base.getindex{PWS1<:PiecewiseSpace,PWS2<:PiecewiseSpace}(d::TensorSpace{@compat(Tuple{PWS1,PWS2})},i::Integer,j::Integer)=d[1][i]⊗d[2][j]
Base.getindex{PWS1<:PiecewiseSpace,PWS2<:PiecewiseSpace}(d::TensorSpace{@compat(Tuple{PWS1,PWS2})},i::Range,j::Range)=PiecewiseSpace(d[1][i])⊗PiecewiseSpace(d[2][j])

## ProductFun

##  Piecewise

function pieces{PS<:PiecewiseSpace}(U::ProductFun{PS})
    ps=space(U,1)
    sp2=space(U,2)
    m=length(ps)
    C=coefficients(U)
    [ProductFun(C[k:m:end,:],ps[k],sp2) for k=1:m]
end
