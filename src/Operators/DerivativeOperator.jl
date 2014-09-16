export DerivativeOperator


immutable DerivativeOperator{T<:Number,S<:FunctionSpace} <: BandedOperator{T}
    space::S        # the domain space
    order::Int
end
DerivativeOperator{S<:PeriodicDomainSpace}(sp::S,k::Integer)=DerivativeOperator{Complex{Float64},S}(sp,k)
DerivativeOperator{S<:FunctionSpace}(sp::S,k::Integer)=DerivativeOperator{Float64,S}(sp,k)

DerivativeOperator(sp::FunctionSpace)=DerivativeOperator(sp,1)
DerivativeOperator()=DerivativeOperator(AnySpace())
DerivativeOperator(k::Integer)=DerivativeOperator(AnySpace(),k)

DerivativeOperator(d::IntervalDomain)=DerivativeOperator(ChebyshevSpace(d))
DerivativeOperator(d::PeriodicDomain)=DerivativeOperator(LaurentSpace(d))

domain(D::DerivativeOperator)=domain(D.space)

domainspace(M::DerivativeOperator)=M.space


## Overrideable
bandinds{T,S<:IntervalDomainSpace}(D::DerivativeOperator{T,S})=0,D.order
bandinds{T,S<:PeriodicDomainSpace}(D::DerivativeOperator{T,S})=0,0





#promoting domain space is allowed to change range space
promotedomainspace(D::DerivativeOperator,sp::AnySpace)=D
function promotedomainspace{S<:FunctionSpace}(D::DerivativeOperator,sp::S)
    if domain(sp) == AnyDomain()
         DerivativeOperator(S(domain(D)),D.order)
    else
        DerivativeOperator(sp,D.order)
    end
end


## simplify higher order derivatives
function *(D1::DerivativeOperator,D2::DerivativeOperator)
    @assert domain(D1) == domain(D2)
    
    DerivativeOperator(D2.space,D1.order+D2.order)
end

#^(D1::DerivativeOperator,k::Integer)=DerivativeOperator(D1.order*k,D1.space)



