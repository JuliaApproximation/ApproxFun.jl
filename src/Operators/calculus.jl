export DerivativeOperator,IntegrationOperator


for TT in (:DerivativeOperator,:IntegrationOperator)
    @eval begin
        immutable $TT{T<:Number,S<:FunctionSpace} <: BandedOperator{T}
            space::S        # the domain space
            order::Int
        end
        $TT{S<:PeriodicDomainSpace}(sp::S,k::Integer)=$TT{Complex{Float64},S}(sp,k)
        $TT{S<:FunctionSpace}(sp::S,k::Integer)=$TT{Float64,S}(sp,k)
        
        $TT(sp::FunctionSpace)=$TT(sp,1)
        $TT()=$TT(AnySpace())
        $TT(k::Integer)=$TT(AnySpace(),k)
        
        $TT(d::PeriodicDomain,n::Integer)=$TT(LaurentSpace(d),n)
        $TT(d::Domain)=$TT(d,1)

        
        domain(D::$TT)=domain(D.space)
        
        domainspace(M::$TT)=M.space
        bandinds{T,S<:PeriodicDomainSpace}(D::$TT{T,S})=0,0
    end
end

# the default domain space is higher to avoid negative ultraspherical spaces
DerivativeOperator(d::IntervalDomain,n::Integer)=DerivativeOperator(ChebyshevSpace(d),n)
IntegrationOperator(d::IntervalDomain,n::Integer)=IntegrationOperator(UltrasphericalSpace{1}(d),n)


#promoting domain space is allowed to change range space
# for integration, we fall back on existing conversion for now
promotedomainspace(D::DerivativeOperator,sp::AnySpace)=D

function promotedomainspace{S<:FunctionSpace}(D::DerivativeOperator,sp::S)
    if domain(sp) == AnyDomain()
         DerivativeOperator(S(domain(D)),D.order)
    else
        DerivativeOperator(sp,D.order)
    end
end        


## simplify higher order derivatives/integration
function *(D1::DerivativeOperator,D2::DerivativeOperator)
    @assert domain(D1) == domain(D2)
    
    DerivativeOperator(D2.space,D1.order+D2.order)
end


## Overrideable
bandinds{T,S<:IntervalDomainSpace}(D::DerivativeOperator{T,S})=0,D.order
bandinds{T,S<:IntervalDomainSpace}(D::IntegrationOperator{T,S})=-D.order,0





## Convenience routines

Base.diff(d::DomainSpace,μ::Integer)=DerivativeOperator(d,μ)
Base.diff(d::Domain,μ::Integer)=DerivativeOperator(d,μ)
Base.diff(d::Domain)=Base.diff(d,1)

integrate(d::Domain)=IntegrationOperator(d,1)



#^(D1::DerivativeOperator,k::Integer)=DerivativeOperator(D1.order*k,D1.space)



