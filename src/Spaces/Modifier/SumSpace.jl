export ⊕

## SumSpace{T,S,V} encodes a space that can be decoupled as f(x) = a(x) + b(x) where a is in S and b is in V


immutable SumSpace{S<:FunctionSpace,V<:FunctionSpace,T<:Number,D<:Domain} <: FunctionSpace{T,D}
    spaces::(S,V)
    SumSpace(d::Domain)=new((S(d),V(d)))
    SumSpace(sp::(S,V))=new(sp)
end

function SumSpace{T<:Number,D}(A::(FunctionSpace{T,D},FunctionSpace{T,D}))
    @assert domain(A[1])==domain(A[2])
    SumSpace{typeof(A[1]),typeof(A[2]),T,D}(A)
end

SumSpace(A::FunctionSpace,B::FunctionSpace)=SumSpace((A,B))


typealias PeriodicSumSpace{S,V,T} SumSpace{S,V,T,PeriodicInterval}
typealias IntervalSumSpace{S,V,T} SumSpace{S,V,T,Interval}




⊕(A::FunctionSpace,B::FunctionSpace)=SumSpace(A,B)
⊕(f::Fun,g::Fun)=Fun(interlace(coefficients(f),coefficients(g)),space(f)⊕space(g))



Base.getindex(S::SumSpace,k)=S.spaces[k]

domain(A::SumSpace)=domain(A[1])
evaluate{D<:SumSpace,T}(f::Fun{D,T},x)=evaluate(Fun(f.coefficients[1:2:end],space(f)[1]),x)+evaluate(Fun(f.coefficients[2:2:end],space(f)[2]),x)


spacescompatible{S,T}(A::SumSpace{S,T},B::SumSpace{S,T})=spacescompatible(A.spaces[1],B[1]) && spacescompatible(A.spaces[2],B[2])





## calculus

# assume first domain has 1 as a basis element

Base.ones{T<:Number}(::Type{T},S::SumSpace)=ones(T,S[1])⊕zeros(T,S[2])
Base.ones(S::SumSpace)=ones(S[1])⊕zeros(S[2])


# vec

Base.vec(S::SumSpace)=S.spaces
Base.vec{S<:SumSpace,T}(f::Fun{S,T})=Fun[Fun(f.coefficients[j:2:end],space(f)[j]) for j=1:2]



## values

itransform(S::SumSpace,cfs)=Fun(cfs,S)[points(S,length(cfs))]


