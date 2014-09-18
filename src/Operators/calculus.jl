export Derivative,Integral


for TT in (:Derivative,:Integral)
    @eval begin
        immutable $TT{S<:FunctionSpace,T<:Number} <: BandedOperator{T}
            space::S        # the domain space
            order::Int
        end
        $TT{S<:PeriodicDomainSpace}(sp::S,k::Integer)=$TT{S,Complex{Float64}}(sp,k)
        $TT{S<:FunctionSpace}(sp::S,k::Integer)=$TT{S,Float64}(sp,k)
        
        $TT(sp::FunctionSpace)=$TT(sp,1)
        $TT()=$TT(AnySpace())
        $TT(k::Integer)=$TT(AnySpace(),k)
        
        $TT(d::PeriodicDomain,n::Integer)=$TT(LaurentSpace(d),n)
        $TT(d::Domain)=$TT(d,1)

        
        domain(D::$TT)=domain(D.space)
        
        domainspace(D::$TT)=D.space
        rangespace{T,S<:PeriodicDomainSpace}(D::$TT{S,T})=D.space        #assume rangespace is the same
        bandinds{T,S<:PeriodicDomainSpace}(D::$TT{S,T})=0,0
    end
end

# the default domain space is higher to avoid negative ultraspherical spaces
Derivative(d::IntervalDomain,n::Integer)=Derivative(ChebyshevSpace(d),n)
Integral(d::IntervalDomain,n::Integer)=Integral(UltrasphericalSpace{1}(d),n)


#promoting domain space is allowed to change range space
# for integration, we fall back on existing conversion for now
promotedomainspace(D::Derivative,sp::AnySpace)=D

function promotedomainspace{S<:FunctionSpace}(D::Derivative,sp::S)
    if domain(sp) == AnyDomain()
         Derivative(S(domain(D)),D.order)
    else
        Derivative(sp,D.order)
    end
end        


## simplify higher order derivatives/integration
function *(D1::Derivative,D2::Derivative)
    @assert domain(D1) == domain(D2)
    
    Derivative(D2.space,D1.order+D2.order)
end


## Overrideable
bandinds{T,S<:IntervalDomainSpace}(D::Derivative{S,T})=0,D.order
bandinds{T,S<:IntervalDomainSpace}(D::Integral{S,T})=-D.order,0





## Convenience routines

Base.diff(d::DomainSpace,μ::Integer)=Derivative(d,μ)
Base.diff(d::Domain,μ::Integer)=Derivative(d,μ)
Base.diff(d::Domain)=Base.diff(d,1)

integrate(d::Domain)=Integral(d,1)



#^(D1::Derivative,k::Integer)=Derivative(D1.order*k,D1.space)



