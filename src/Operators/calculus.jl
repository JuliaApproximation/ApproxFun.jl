export Derivative,Integral


abstract CalculusOperator{S,OT,T}<:BandedOperator{T}


## Note that all functions called in calculus_operator must be exported



macro calculus_operator(Op)
    AbstOp=parse("Abstract"*string(Op))
    WrappOp=parse(string(Op)*"Wrapper")
    return esc(quote
        # The SSS, TTT are to work around #9312
        abstract $AbstOp{SSS,OT,TTT} <: CalculusOperator{SSS,OT,TTT}

        immutable $Op{S<:FunctionSpace,OT,T} <: $AbstOp{S,OT,T}
            space::S        # the domain space
            order::OT
        end
        immutable $WrappOp{BT<:BandedOperator,S<:FunctionSpace,OT,T} <: $AbstOp{S,OT,T}
            op::BT
            order::OT
        end


        ## Constructors
        $Op{T}(::Type{T},sp::FunctionSpace,k)=$Op{typeof(sp),typeof(k),T}(sp,k)
        $Op(::Type{Any},sp::FunctionSpace,k)=$Op(sp,k)


        $Op(sp::FunctionSpace{ComplexBasis},k)=$Op{typeof(sp),typeof(k),Complex{real(eltype(domain(sp)))}}(sp,k)
        $Op(sp::FunctionSpace,k)=$Op{typeof(sp),typeof(k),eltype(domain(sp))}(sp,k)

        $Op(sp::FunctionSpace)=$Op(sp,1)
        $Op()=$Op(UnsetSpace())
        $Op(k::Integer)=$Op(UnsetSpace(),k)

        $Op(d::Domain,n)=$Op(Space(d),n)
        $Op(d::Domain)=$Op(d,1)
        $Op(d::Vector)=$Op(Space(d),1)
        $Op(d::Vector,n)=$Op(Space(d),n)

        Base.convert{BT<:$Op}(::Type{BT},D::BT)=D
        function Base.convert{BT<:Operator}(::Type{BT},D::$Op)
            if eltype(BT)==eltype(D)
                D
            else
                $Op{typeof(D.space),typeof(D.order),eltype(BT)}(D.space,D.order)
            end
        end

        $WrappOp{T<:Number}(op::BandedOperator{T},order)=$WrappOp{typeof(op),typeof(domainspace(op)),typeof(order),T}(op,order)
        $WrappOp{T<:Number}(op::BandedOperator{T})=$WrappOp(op,1)

        Base.convert{BT<:$WrappOp}(::Type{BT},D::BT)=D
        function Base.convert{BT<:Operator}(::Type{BT},D::$WrappOp)
            if eltype(BT)==eltype(D)
                D
            else
                $WrappOp(convert(BandedOperator{eltype(BT)},D.op),D.order)
            end
        end

        ## Routines
        domain(D::$Op)=domain(D.space)
        domainspace(D::$Op)=D.space

        addentries!{OT,T}(::$Op{UnsetSpace,OT,T},A,kr::Range)=error("Spaces cannot be inferred for operator")

        function addentries!{S,OT,T}(D::$Op{S,OT,T},A,kr::Range)
            # Default is to convert to Canonical and apply operator there
            sp=domainspace(D)
            csp=canonicalspace(sp)
            if conversion_type(csp,sp)==csp   # Conversion(sp,csp) is not banded, or sp==csp
                error("Override "*string($Op)*"(::"*string(typeof(sp))*","*string(D.order)*")")
            end
            addentries!(TimesOperator([$Op(csp,D.order),Conversion(sp,csp)]),A,kr)
        end

        function bandinds(D::$Op)
            sp=domainspace(D)
            csp=canonicalspace(sp)
            if conversion_type(csp,sp)==csp   # Conversion(sp,csp) is not banded, or sp==csp
                error("Override "*string($Op)*"(::"*string(typeof(sp))*","*string(D.order)*")")
            end
            bandinds(TimesOperator([$Op(csp,D.order),Conversion(sp,csp)]))
        end

        # corresponds to default implementation
        function rangespace{S,T}(D::$Op{S,T})
            sp=domainspace(D)
            csp=canonicalspace(sp)
            if conversion_type(csp,sp)==csp   # Conversion(sp,csp) is not banded, or sp==csp
                error("Override *"string($Op)*"(::"*string(typeof(sp))*","*string(D.order)*")")
            end
            rangespace($Op(canonicalspace(domainspace(D)),D.order))
        end
        rangespace{T}(D::$Op{UnsetSpace,T})=UnsetSpace()

        #promoting domain space is allowed to change range space
        # for integration, we fall back on existing conversion for now
        promotedomainspace(D::$AbstOp,sp::UnsetSpace)=D
        promotedomainspace(D::$AbstOp,sp::AnySpace)=D


        function promotedomainspace{S<:FunctionSpace}(D::$AbstOp,sp::S)
            if isambiguous(domain(sp))
                $Op(S(domain(D)),D.order)
            else
                $Op(sp,D.order)
            end
        end

        choosedomainspace(M::$Op{UnsetSpace},sp)=sp  # we assume


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



@calculus_operator(Derivative)
@calculus_operator(Integral)

for (ATYP,TYP) in ((:AbstractDerivative,:Derivative),(:AbstractIntegral,:Integral))
    @eval begin
        function *(D1::$ATYP,D2::$ATYP)
            @assert domain(D1) == domain(D2)

            $TYP(domainspace(D2),D1.order+D2.order)
        end
    end
end






## Overrideable


## Convenience routines


integrate(d::Domain)=Integral(d,1)


# Default is to use ops
differentiate{S,T}(f::Fun{S,T})=Derivative(space(f))*f
integrate{S,T}(f::Fun{S,T})=Integral(space(f))*f


#^(D1::Derivative,k::Integer)=Derivative(D1.order*k,D1.space)






# Multivariate



@calculus_operator(Laplacian)

Laplacian(S::FunctionSpace,k)=Laplacian{typeof(S),Int,BandedMatrix{eltype(S)}}(S,k)
Laplacian(S)=Laplacian(S,1)
