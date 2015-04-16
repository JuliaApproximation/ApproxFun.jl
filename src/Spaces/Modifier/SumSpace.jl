export ⊕

## SumSpace{T,S,V} encodes a space that can be decoupled as f(x) = a(x) + b(x) where a is in S and b is in V


immutable SumSpace{S<:FunctionSpace,V<:FunctionSpace,T,d} <: FunctionSpace{T,d}
    spaces::(S,V)
    SumSpace(dom::Domain)=new((S(dom),V(dom)))
    SumSpace(sp::(S,V))=new(sp)
end

SumSpace{T1,T2,d}(A::(FunctionSpace{T1,d},FunctionSpace{T2,d}))=SumSpace{typeof(A[1]),typeof(A[2]),promote_type(T1,T2),d}(A)

SumSpace(A::FunctionSpace,B::FunctionSpace)=SumSpace((A,B))






⊕(A::FunctionSpace,B::FunctionSpace)=SumSpace(A,B)
⊕(f::Fun,g::Fun)=Fun(interlace(coefficients(f),coefficients(g)),space(f)⊕space(g))




Base.getindex(S::SumSpace,k)=S.spaces[k]

domain(A::SumSpace)=domain(A[1])



spacescompatible{S,T,d}(A::SumSpace{S,T,d},B::SumSpace{S,T,d})=spacescompatible(A.spaces[1],B[1]) && spacescompatible(A.spaces[2],B[2])





## routines

evaluate{D<:SumSpace,T}(f::Fun{D,T},x)=evaluate(vec(f,1),x)+evaluate(vec(f,2),x)
for OP in (:differentiate,:integrate)
    @eval $OP{D<:SumSpace,T}(f::Fun{D,T})=$OP(vec(f,1))⊕$OP(vec(f,2))
end

# assume first domain has 1 as a basis element

Base.ones{T<:Number}(::Type{T},S::SumSpace)=ones(T,S[1])⊕zeros(T,S[2])
Base.ones(S::SumSpace)=ones(S[1])⊕zeros(S[2])


# vec

Base.vec{D<:SumSpace,T}(f::Fun{D,T},k)=k==1?Fun(f.coefficients[1:2:end],space(f)[1]):Fun(f.coefficients[2:2:end],space(f)[2])
Base.vec(S::SumSpace)=S.spaces
Base.vec{S<:SumSpace,T}(f::Fun{S,T})=Fun[vec(f,j) for j=1:2]



## values

itransform(S::SumSpace,cfs)=Fun(cfs,S)[points(S,length(cfs))]


