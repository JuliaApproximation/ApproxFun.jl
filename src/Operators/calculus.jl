export Derivative,Integral

abstract AbstractDerivative{T} <: BandedOperator{T}
abstract AbstractIntegral{T} <: BandedOperator{T}

immutable Derivative{S<:FunctionSpace,T<:Number} <: AbstractDerivative{T}
    space::S        # the domain space
    order::Int
end

immutable Integral{S<:FunctionSpace,T<:Number} <: AbstractIntegral{T}
    space::S        # the domain space
    order::Int
end

for TT in (:Derivative,:Integral)
    @eval begin
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
        
        addentries!{T}(::$TT{AnySpace,T},A::ShiftArray,kr::Range)=error("Spaces cannot be inferred for operator")
        
        function addentries!{S,T}(D::$TT{S,T},A::ShiftArray,kr::Range)   
            # Default is to convert to Canonical and d
            sp=domainspace(D)
            csp=canonicalspace(sp)
            addentries!(TimesOperator([$TT(csp,D.order),Conversion(sp,csp)]),A,kr)
        end
        
        function bandinds(D::$TT)
            sp=domainspace(D)
            csp=canonicalspace(sp)
            bandinds(TimesOperator([$TT(csp,D.order),Conversion(sp,csp)])) 
        end
        
        rangespace{S,T}(D::$TT{S,T})=rangespace($TT(canonicalspace(domainspace(D)),D.order))
        rangespace{T}(D::$TT{AnySpace,T})=AnySpace()
    end
end

# the default domain space is higher to avoid negative ultraspherical spaces
Derivative(d::IntervalDomain,n::Integer)=Derivative(ChebyshevSpace(d),n)
Integral(d::IntervalDomain,n::Integer)=Integral(UltrasphericalSpace{1}(d),n)


#promoting domain space is allowed to change range space
# for integration, we fall back on existing conversion for now
promotedomainspace(D::AbstractDerivative,sp::AnySpace)=D

function promotedomainspace{S<:FunctionSpace}(D::AbstractDerivative,sp::S)
    if domain(sp) == AnyDomain()
         Derivative(S(domain(D)),D.order)
    else
        Derivative(sp,D.order)
    end
end        


## simplify higher order derivatives/integration
function *(D1::AbstractDerivative,D2::AbstractDerivative)
    @assert domain(D1) == domain(D2)
    
    Derivative(domainspace(D2),D1.order+D2.order)
end


## Overrideable






## Convenience routines

Base.diff(d::DomainSpace,μ::Integer)=Derivative(d,μ)
Base.diff(d::Domain,μ::Integer)=Derivative(d,μ)
Base.diff(d::Domain)=Base.diff(d,1)

integrate(d::Domain)=Integral(d,1)


# Default is to use ops
differentiate(f::Fun)=Derivative(space(f))*f


#^(D1::Derivative,k::Integer)=Derivative(D1.order*k,D1.space)



