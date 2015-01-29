export ConstantOperator, BasisFunctional


immutable ConstantOperator{T<:Number} <: BandedOperator{T}
    c::T
end

ConstantOperator(L::UniformScaling)=ConstantOperator(L.λ)
IdentityOperator()=ConstantOperator(1.0)

bandinds(T::ConstantOperator)=0,0

addentries!(C::ConstantOperator,A,kr::Range)=toeplitz_addentries!([.5C.c],A,kr)

==(C1::ConstantOperator,C2::ConstantOperator)=C1.c==C2.c


Base.convert{T<:Number}(::Type{BandedOperator{T}},C::ConstantOperator)=ConstantOperator{T}(C.c)

## Algebra

for op in (:+,:-,:*)
    @eval ($op)(A::ConstantOperator,B::ConstantOperator)=ConstantOperator($op(A.c,B.c))
end


## Basis Functional

immutable BasisFunctional <: Functional{Float64}
    k::Integer
end


Base.getindex(op::BasisFunctional,k::Integer)=(k==op.k)?1.:0.
Base.getindex(op::BasisFunctional,k::Range1)=convert(Vector{Float64},k.==op.k)

immutable FillFunctional{T<:Number} <: Functional{T}
    c::T
end


Base.getindex(op::FillFunctional,k::Integer)=op.c
Base.getindex(op::FillFunctional,k::Range)=fill(op.c,length(k))

## Zero is a special operator: it makes sense on all spaces, and between all spaces

immutable ZeroOperator{T,S,V} <: BandedOperator{T}
    domainspace::S
    rangespace::V
end

ZeroOperator{T<:Number,S,V}(::Type{T},d::S,v::V)=ZeroOperator{T,S,V}(d,v)
ZeroOperator{S,V}(d::S,v::V)=ZeroOperator(Float64,d,v)
ZeroOperator()=ZeroOperator(AnySpace(),AnySpace())
ZeroOperator{T<:Number}(::Type{T})=ZeroOperator(T,AnySpace(),AnySpace())

domainspace(Z::ZeroOperator)=Z.domainspace
rangespace(Z::ZeroOperator)=Z.rangespace

bandinds(T::ZeroOperator)=0,0

addentries!(C::ZeroOperator,A,kr::Range)=A

promotedomainspace(Z::ZeroOperator,sp::AnySpace)=Z
promoterangespace(Z::ZeroOperator,sp::AnySpace)=Z
promotedomainspace(Z::ZeroOperator,sp::UnsetSpace)=Z
promoterangespace(Z::ZeroOperator,sp::UnsetSpace)=Z
promotedomainspace(Z::ZeroOperator,sp::FunctionSpace)=ZeroOperator(sp,rangespace(Z))
promoterangespace(Z::ZeroOperator,sp::FunctionSpace)=ZeroOperator(domainspace(Z),sp)
promotedomainspace(Z::ZeroOperator,sp::FunctionSpace)=ZeroOperator(sp,rangespace(Z))
promoterangespace(Z::ZeroOperator,sp::FunctionSpace)=ZeroOperator(domainspace(Z),sp)


immutable ZeroFunctional{S<:FunctionSpace,T<:Number} <: Functional{T}
    domainspace::S
end
ZeroFunctional(sp::FunctionSpace)=ZeroFunctional{typeof(sp),Float64}(sp)
ZeroFunctional{T<:Number}(::Type{T},sp::FunctionSpace)=ZeroFunctional{typeof(sp),T}(sp)
ZeroFunctional{T<:Number}(::Type{T})=ZeroFunctional(T,AnySpace())
ZeroFunctional()=ZeroFunctional(AnySpace())

Base.convert{T}(::Type{Functional{T}},Z::ZeroFunctional)=ZeroFunctional(T,Z.domainspace)

domainspace(Z::ZeroFunctional)=Z.domainspace
promotedomainspace(Z::ZeroFunctional,sp::FunctionSpace)=ZeroFunctional(sp)

Base.getindex{T}(op::ZeroFunctional{T},k::Integer)=zero(T)
Base.getindex{T}(op::ZeroFunctional{T},k::Range1)=zeros(T,length(k))



