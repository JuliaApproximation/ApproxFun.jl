export Derivative,Integral


abstract CalculusOperator{S,T}<:BandedOperator{T}

macro calculus_operator(Op,AbstOp,WrappOp)
    return esc(quote        
        # The SSS, TTT are to work around #9312
        abstract $AbstOp{SSS,TTT} <: CalculusOperator{SSS,TTT} 
    
        immutable $Op{S<:FunctionSpace,T<:Number} <: $AbstOp{S,T}
            space::S        # the domain space
            order::Int
        end                       
        immutable $WrappOp{BT<:BandedOperator,S<:FunctionSpace,T<:Number} <: $AbstOp{S,T}
            op::BT
            order::Int
        end    

            
        ## Constructors  
        $Op{T}(::Type{T},sp::FunctionSpace,k)=$Op{typeof(sp),T}(sp,k)
              
        $Op(sp::FunctionSpace{RealBasis},k)=$Op{typeof(sp),Float64}(sp,k)
        $Op(sp::FunctionSpace{ComplexBasis},k)=$Op{typeof(sp),Complex{Float64}}(sp,k)        
        
        $Op(sp::FunctionSpace)=$Op(sp,1)
        $Op()=$Op(AnySpace())
        $Op(k::Integer)=$Op(AnySpace(),k)
        
        $Op(d::Domain,n)=$Op(Space(d),n)
        $Op(d::Domain)=$Op(d,1)
        $Op(d::Vector)=$Op(Space(d),1)
        $Op(d::Vector,n)=$Op(Space(d),n)        
        
        Base.convert{T}(::Type{BandedOperator{T}},D::$Op)=$Op(T,D.space,D.order)
        
        $WrappOp{T<:Number}(op::BandedOperator{T},order::Integer)=$WrappOp{typeof(op),typeof(domainspace(op)),T}(op,order)
        $WrappOp{T<:Number}(op::BandedOperator{T})=$WrappOp(op,1)        
        
        ## Routines
        domain(D::$Op)=domain(D.space)       
        domainspace(D::$Op)=D.space
        
        addentries!{T}(::$Op{AnySpace,T},A,kr::Range)=error("Spaces cannot be inferred for operator")
        
        function addentries!{S,T}(D::$Op{S,T},A,kr::Range)   
            # Default is to convert to Canonical and d
            sp=domainspace(D)
            csp=canonicalspace(sp)
            if csp==sp
                error("Override "*string($Op)*"(::"*string(typeof(sp))*","*string(D.order)*")")
            end
            addentries!(TimesOperator([$Op(csp,D.order),Conversion(sp,csp)]),A,kr)
        end
        
        function bandinds(D::$Op)
            sp=domainspace(D)
            csp=canonicalspace(sp)
            if csp==sp
                error("Override "*string($Op)*"(::"*string(typeof(sp))*","*string(D.order)*")")
            end     
            bandinds(TimesOperator([$Op(csp,D.order),Conversion(sp,csp)])) 
        end

        # corresponds to default implementation        
        function rangespace{S,T}(D::$Op{S,T})
            sp=domainspace(D)
            csp=canonicalspace(sp)
            if csp==sp
                error("Override *"string($Op)*"(::"*string(typeof(sp))*","*string(D.order)*")")
            end      
            rangespace($Op(canonicalspace(domainspace(D)),D.order))
        end
        rangespace{T}(D::$Op{AnySpace,T})=AnySpace()     
        
        #promoting domain space is allowed to change range space
        # for integration, we fall back on existing conversion for now
        promotedomainspace(D::$AbstOp,sp::AnySpace)=D
        
        function promotedomainspace{S<:FunctionSpace}(D::$AbstOp,sp::S)
            if domain(sp) == AnyDomain()
                $Op(S(domain(D)),D.order)
            else
                $Op(sp,D.order)
            end
        end
        
        #Wrapper just adds the operator it wraps
        addentries!(D::$WrappOp,A,k::Range)=addentries!(D.op,A,k)          
        rangespace(D::$WrappOp)=rangespace(D.op)
        domainspace(D::$WrappOp)=domainspace(D.op)        
        bandinds(D::$WrappOp)=bandinds(D.op)        
    end)
#     for func in (:rangespace,:domainspace,:bandinds)
#         # We assume the operator wrapped has the correct spaces
#         @eval $func(D::$WrappOp)=$func(D.op)
#     end 
end



@calculus_operator(Derivative,AbstractDerivative,DerivativeWrapper)
@calculus_operator(Integral,AbstractIntegral,IntegralWrapper)

      



## simplify higher order derivatives/integration
function *(D1::AbstractDerivative,D2::AbstractDerivative)
    @assert domain(D1) == domain(D2)
    
    Derivative(domainspace(D2),D1.order+D2.order)
end


## Overrideable


## Convenience routines

Base.diff(d::FunctionSpace,μ::Integer)=Derivative(d,μ)
Base.diff(d::Domain,μ::Integer)=Derivative(d,μ)
Base.diff(d::Domain)=Base.diff(d,1)

integrate(d::Domain)=Integral(d,1)


# Default is to use ops
differentiate{S,T}(f::Fun{S,T})=Derivative(space(f))*f
integrate{S,T}(f::Fun{S,T})=Integral(space(f))*f


#^(D1::Derivative,k::Integer)=Derivative(D1.order*k,D1.space)






