export ConstantOperator, IdentityOperator, BasisFunctional

struct ConstantOperator{T,DS} <: Operator{T}
    λ::T
    space::DS
    (::Type{ConstantOperator{T,DS}}){T,DS}(c::Number,sp::DS) = new{T,DS}(convert(T,c),sp)
    (::Type{ConstantOperator{T,DS}}){T,DS}(L::UniformScaling,sp::DS) = new{T,DS}(convert(T,L.λ),sp)
end


ConstantOperator{T}(::Type{T},c,sp::Space) = ConstantOperator{T,typeof(sp)}(c,sp)
ConstantOperator{T}(::Type{T},c) = ConstantOperator(T,c,UnsetSpace())
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

bandinds(T::ConstantOperator) = 0,0

isbandedblockbanded(::ConstantOperator) = true
blockbandinds(::ConstantOperator) = 0,0
subblockbandinds(::ConstantOperator,k::Integer) = 0

getindex(C::ConstantOperator,k::Integer,j::Integer) =
    k==j?eltype(C)(C.λ):zero(eltype(C))


==(C1::ConstantOperator,C2::ConstantOperator) = C1.λ==C2.λ

function convert{T}(::Type{Operator{T}},C::ConstantOperator)
    if T == eltype(C)
        C
    else
        ConstantOperator{T,typeof(C.space)}(C.λ,C.space)
    end
end

# zero needs to be different since it can take a space to
# a ConstantSpace, in creating functionals
convert{T}(::Type{Operator{T}},x::Number) =
    x==0 ? ZeroOperator(T) : Multiplication(T(x))
convert{T}(::Type{Operator{T}},L::UniformScaling) =
    ConstantOperator(T,L.λ)

convert(::Type{Operator},n::Number) = Operator{typeof(n)}(n)
convert(::Type{Operator},L::UniformScaling) = ConstantOperator(L.λ)

## Algebra

for op in (:+,:-,:*)
    @eval ($op)(A::ConstantOperator,B::ConstantOperator) = ConstantOperator($op(A.λ,B.λ))
end


## Basis Functional

struct BasisFunctional{T} <: Operator{T}
    k::Integer
end

@functional BasisFunctional

BasisFunctional(k) = BasisFunctional{Float64}(k)

bandinds(B::BasisFunctional) = 0,B.k-1
domainspace(B::BasisFunctional) = ℓ⁰

convert{T}(::Type{Operator{T}},B::BasisFunctional) = BasisFunctional{T}(B.k)

Base.getindex{T}(op::BasisFunctional{T},k::Integer) = (k==op.k)?one(T):zero(T)
Base.getindex{T}(op::BasisFunctional{T},k::Range) = convert(Vector{T},k.==op.k)

struct FillFunctional{T} <: Operator{T}
    λ::T
end

@functional FillFunctional

domainspace(B::FillFunctional) = ℓ⁰

Base.getindex(op::FillFunctional,k::Integer)=op.λ
Base.getindex(op::FillFunctional,k::Range)=fill(op.λ,length(k))

## Zero is a special operator: it makes sense on all spaces, and between all spaces

struct ZeroOperator{T,S,V} <: Operator{T}
    domainspace::S
    rangespace::V
end

ZeroOperator{T}(::Type{T},d::Space,v::Space) = ZeroOperator{T,typeof(d),typeof(v)}(d,v)
ZeroOperator{T}(::Type{T},S::Space) = ZeroOperator(T,S,ZeroSpace(S))
ZeroOperator(d::Space,v::Space) = ZeroOperator(Float64,d,v)
ZeroOperator(S::Space) = ZeroOperator(S,ZeroSpace(S))
ZeroOperator() = ZeroOperator(UnsetSpace(),UnsetSpace())
ZeroOperator{T}(::Type{T}) = ZeroOperator(T,UnsetSpace(),UnsetSpace())


convert{T}(::Type{Operator{T}},Z::ZeroOperator) =
    ZeroOperator(T,Z.domainspace,Z.rangespace)


Base.zero(A::Operator) = ZeroOperator(domainspace(A),rangespace(A))


domainspace(Z::ZeroOperator)=Z.domainspace
rangespace(Z::ZeroOperator)=Z.rangespace

bandinds(T::ZeroOperator) = 720,-720   # 6!
blockbandinds(T::ZeroOperator) = 720,-720   # 6!
subblockbandinds(T::ZeroOperator,k::Integer) = k == 1 ? 720 : -720

isbandedblockbandedabove(::ZeroOperator) = true
isbandedblockbandedbelow(::ZeroOperator) = true

getindex(C::ZeroOperator,k::Integer) = zero(eltype(C))
getindex(C::ZeroOperator,k::Integer,j::Integer) = zero(eltype(C))

promotedomainspace(Z::ZeroOperator,sp::UnsetSpace) = Z
promoterangespace(Z::ZeroOperator,sp::UnsetSpace) = Z
promotedomainspace(Z::ZeroOperator,sp::Space) = ZeroOperator(sp,rangespace(Z))
promoterangespace(Z::ZeroOperator,sp::Space) = ZeroOperator(domainspace(Z),sp)




isconstop(::Union{ZeroOperator,ConstantOperator})=true
isconstop(S::SpaceOperator)=isconstop(S.op)
isconstop(::)=false

iszeroop(::ZeroOperator) = true
iszeroop(A::ConstantOperator) = A.λ==0.0
iszeroop(A) = false

convert{T<:Number}(::Type{T},::ZeroOperator) = zero(T)
convert{T<:Number}(::Type{T},C::ConstantOperator) = convert(T,C.λ)
convert{T<:Number}(::Type{T},S::SpaceOperator) = convert(T,S.op)



## SubOperator convert
## Special case for ZeroOperator
for (TYP,ZER) in ((:Matrix,:zeros),(:BandedMatrix,:bzeros),(:RaggedMatrix,:rzeros),
                    (:BlockBandedMatrix,:bbzeros))
    @eval convert{T,ZO<:ZeroOperator}(::Type{$TYP},S::SubOperator{T,ZO}) = $ZER(S)
end
