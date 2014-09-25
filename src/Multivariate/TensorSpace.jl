
abstract MultivariateDomain
##TODO: MultivariateDomain{2}
abstract BivariateDomain <: MultivariateDomain





immutable ProductDomain{D<:Domain} <:BivariateDomain
    domains::Vector{D} 
end

ProductDomain(A,B)=ProductDomain([A,B])
*(A::Domain,B::Domain)=ProductDomain(A,B)

Base.transpose(d::ProductDomain)=ProductDomain(d[2],d[1])
Base.getindex(d::ProductDomain,k::Integer)=d.domains[k]


abstract MultivariateFunctionSpace
abstract BivariateFunctionSpace <: MultivariateFunctionSpace


fromcanonical(d::MultivariateDomain,x::Tuple)=fromcanonical(d,x...)
tocanonical(d::MultivariateDomain,x::Tuple)=tocanonical(d,x...)
fromcanonical(d::MultivariateFunctionSpace,x...)=fromcanonical(domain(d),x...)
tocanonical(d::MultivariateFunctionSpace,x...)=tocanonical(domain(d),x...)


# This means x are represented as space S and y are represented as space T
abstract AbstractProductSpace{S,T} <: BivariateFunctionSpace

immutable TensorSpace{S<:FunctionSpace,T<:FunctionSpace} <:AbstractProductSpace{S,T}
    spaces::(S,T)
end

TensorSpace(A,B)=TensorSpace((A,B))
TensorSpace(A::ProductDomain)=TensorSpace(Space(A[1]),Space(A[2]))
⊗(A::FunctionSpace,B::FunctionSpace)=TensorSpace(A,B)
domain(f::TensorSpace)=domain(f.spaces[1])*domain(f.spaces[2])
Space(sp::ProductDomain)=TensorSpace(sp)

Base.getindex(d::TensorSpace,k::Integer)=d.spaces[k]


immutable ProductSpace{S<:FunctionSpace,T<:FunctionSpace} <: AbstractProductSpace{S,T}
    spacesx::Vector{S}
    spacey::T
end

⊗{S<:FunctionSpace}(A::Vector{S},B::FunctionSpace)=ProductSpace(A,B)
domain(f::ProductSpace)=domain(f.spacesx[1])*domain(f.spacesy)

Base.getindex(d::ProductSpace,k::Integer)=k==1?d.spacesx:d.spacey


space(d::AbstractProductSpace,k::Integer)=d[k]


for TT in (:ProductDomain,:TensorSpace)
    @eval Base.transpose(d::$TT)=$TT(d[2],d[1])
end

