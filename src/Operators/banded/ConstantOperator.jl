export ConstantOperator, IdentityOperator, BasisFunctional

struct ConstantOperator{T,DS} <: Operator{T}
    λ::T
    space::DS
    ConstantOperator{T,DS}(c::Number,sp::DS) where {T,DS} = new{T,DS}(convert(T,c),sp)
    ConstantOperator{T,DS}(L::UniformScaling,sp::DS) where {T,DS} = new{T,DS}(convert(T,L.λ),sp)
end


ConstantOperator(::Type{T},c,sp::Space) where {T} = ConstantOperator{T,typeof(sp)}(c,sp)
ConstantOperator(::Type{T},c) where {T} = ConstantOperator(T,c,UnsetSpace())
ConstantOperator(c::Number,sp::Space) = ConstantOperator(typeof(c),c,sp)
ConstantOperator(c::Number) = ConstantOperator(typeof(c),c)
ConstantOperator(L::UniformScaling) = ConstantOperator(L.λ)
ConstantOperator(L::UniformScaling,sp::Space) = ConstantOperator(L.λ,sp)
IdentityOperator() = ConstantOperator(1.0)
IdentityOperator(S::Space) = ConstantOperator(1.0,S)

for OP in (:domainspace,:rangespace)
    @eval $OP(C::ConstantOperator) = C.space
end

promotedomainspace(C::ConstantOperator,sp::Space) = ConstantOperator(C.λ,sp)
promoterangespace(C::ConstantOperator,sp::Space,cursp::UnsetSpace) = ConstantOperator(C.λ,sp)

bandwidths(T::ConstantOperator) = 0,0

isbandedblockbanded(::ConstantOperator) = true
blockbandwidths(::ConstantOperator) = 0,0
subblockbandwidths(::ConstantOperator) = 0,0

getindex(C::ConstantOperator,k::Integer,j::Integer) =
    k==j ? eltype(C)(C.λ) : zero(eltype(C))


==(C1::ConstantOperator, C2::ConstantOperator) = C1.λ==C2.λ

function convert(::Type{Operator{T}}, C::ConstantOperator) where T
    if T == eltype(C)
        C
    else
        ConstantOperator{T,typeof(C.space)}(C.λ,C.space)
    end
end

# zero needs to be different since it can take a space to
# a ConstantSpace, in creating functionals
convert(::Type{Operator{T}}, x::Number) where {T} =
    x==0 ? ZeroOperator(T) : Multiplication(convert(T,x))
convert(::Type{Operator{T}}, L::UniformScaling) where {T} =
    ConstantOperator(T,L.λ)

convert(::Type{Operator},n::Number) = convert(Operator{typeof(n)}, n)
convert(::Type{Operator},L::UniformScaling) = ConstantOperator(L.λ)

## Algebra

for op in (:+,:-,:*)
    @eval ($op)(A::ConstantOperator,B::ConstantOperator) = ConstantOperator($op(A.λ,B.λ))
end


## Basis Functional

struct BasisFunctional{T,DS} <: Operator{T}
    k::Int
    space::DS
end

@functional BasisFunctional

BasisFunctional(k,sp) = BasisFunctional{Float64,typeof(sp)}(k, sp)
BasisFunctional(k) = BasisFunctional(k, ℓ⁰)

bandwidths(B::BasisFunctional) = 0,B.k-1
domainspace(B::BasisFunctional) = B.space

convert(::Type{Operator{T}},B::BasisFunctional) where {T} = BasisFunctional{T,typeof(B.space)}(B.k,B.space)

Base.getindex(op::BasisFunctional{T},k::Integer) where {T} = (k==op.k) ? one(T) : zero(T)
Base.getindex(op::BasisFunctional{T},k::AbstractRange) where {T} = convert(Vector{T},k.==op.k)

struct FillFunctional{T} <: Operator{T}
    λ::T
end

@functional FillFunctional

domainspace(B::FillFunctional) = ℓ⁰

Base.getindex(op::FillFunctional,k::Integer)=op.λ
Base.getindex(op::FillFunctional,k::AbstractRange)=fill(op.λ,length(k))

## Zero is a special operator: it makes sense on all spaces, and between all spaces

struct ZeroOperator{T,S,V} <: Operator{T}
    domainspace::S
    rangespace::V
end

ZeroOperator(::Type{T},d::Space,v::Space) where {T} = ZeroOperator{T,typeof(d),typeof(v)}(d,v)
ZeroOperator(::Type{T},S::Space) where {T} = ZeroOperator(T,S,ZeroSpace(S))
ZeroOperator(d::Space,v::Space) = ZeroOperator(Float64,d,v)
ZeroOperator(S::Space) = ZeroOperator(S,ZeroSpace(S))
ZeroOperator() = ZeroOperator(UnsetSpace(),UnsetSpace())
ZeroOperator(::Type{T}) where {T} = ZeroOperator(T,UnsetSpace(),UnsetSpace())


convert(::Type{Operator{T}},Z::ZeroOperator) where {T} =
    ZeroOperator(T,Z.domainspace,Z.rangespace)


Base.zero(A::Operator) = ZeroOperator(domainspace(A),rangespace(A))


domainspace(Z::ZeroOperator)=Z.domainspace
rangespace(Z::ZeroOperator)=Z.rangespace

bandwidths(T::ZeroOperator) = -720,-720   # 6!
blockbandwidths(T::ZeroOperator) = -720,-720   # 6!
subblockbandwidths(T::ZeroOperator) =  -720, -720

isbandedblockbandedabove(::ZeroOperator) = true
isbandedblockbandedbelow(::ZeroOperator) = true

getindex(C::ZeroOperator,k::Integer) = zero(eltype(C))
getindex(C::ZeroOperator,k::Integer,j::Integer) = zero(eltype(C))

promotedomainspace(Z::ZeroOperator,sp::UnsetSpace) = Z
promoterangespace(Z::ZeroOperator,sp::UnsetSpace) = Z
promotedomainspace(Z::ZeroOperator,sp::Space) = ZeroOperator(sp,rangespace(Z))
promoterangespace(Z::ZeroOperator,sp::Space) = ZeroOperator(domainspace(Z),sp)




isconstop(::Union{ZeroOperator,ConstantOperator}) = true
isconstop(S::SpaceOperator) = isconstop(S.op)
isconstop(_) = false

iszeroop(::ZeroOperator) = true
iszeroop(A::ConstantOperator) = A.λ==0.0
iszeroop(A) = false

convert(::Type{T},::ZeroOperator) where {T<:Number} = zero(T)
convert(::Type{T},C::ConstantOperator) where {T<:Number} = convert(T,C.λ)
convert(::Type{T},S::SpaceOperator) where {T<:Number} = convert(T,S.op)



## SubOperator convert
## Special case for ZeroOperator
for TYP in (:RaggedMatrix,:Matrix,:BandedMatrix,
            :BlockBandedMatrix,:BandedBlockBandedMatrix)
    @eval $TYP(S::SubOperator{T,ZO}) where {T,ZO<:ZeroOperator} = $TYP(Zeros, S)
end
