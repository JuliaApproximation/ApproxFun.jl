
abstract MultivariateDomain
abstract BivariateDomain <: MultivariateDomain

immutable ProductDomain{D<:Domain} <:BivariateDomain
    domains::Vector{D} 
end

ProductDomain(A,B)=ProductDomain([A,B])
*(A::Domain,B::Domain)=ProductDomain(A,B)


Base.getindex(d::ProductDomain,k::Integer)=d.domains[k]


abstract MultivariateFunctionSpace
abstract BivariateFunctionSpace <: MultivariateFunctionSpace
abstract AbstractProductSpace{S,T} <: BivariateFunctionSpace

immutable TensorSpace{S<:FunctionSpace,T<:FunctionSpace} <:AbstractProductSpace{S,T}
    spaces::(S,T)
end

TensorSpace(A,B)=TensorSpace((A,B))
⊗(A::FunctionSpace,B::FunctionSpace)=TensorSpace(A,B)


Base.getindex(d::TensorSpace,k::Integer)=d.spaces[k]


immutable ProductSpace{S<:FunctionSpace,T<:FunctionSpace} <: AbstractProductSpace{S,T}
    spacesx::Vector{S}
    spacey
end

⊗{S<:FunctionSpace}(A::Vector{S},B::FunctionSpace)=ProductSpace(A,B)

Base.getindex(d::ProductSpace,k::Integer)=k==1?d.spacesx:d.spacey




immutable Disk <: BivariateDomain
    radius::Float64
    center::(Float64,Float64)
end
Disk()=Disk(0.,(0.,0.))

#canonical is rectangle
# we assume radius and centre are zero for now
fromcanonical(D::Disk,x,t)=.5*(1-x)*cos(t),.5*(1-x)*sin(t)
tocanonical(D::Disk,x,y)=1-2sqrt(x^2+y^2),atan2(y,x)


immutable DiskSpace <: BivariateFunctionSpace
    disk::Disk
end

