export Derivative,Integral,Laplacian


abstract CalculusOperator{S,OT,T}<:BandedOperator{T}


## Note that all functions called in calculus_operator must be exported



macro calculus_operator(Op)
    AbstOp=parse("Abstract"*string(Op))
    WrappOp=parse(string(Op)*"Wrapper")
    return esc(quote
        # The SSS, TTT are to work around #9312
        abstract $AbstOp{SSS,OT,TTT} <: CalculusOperator{SSS,OT,TTT}

        immutable $Op{S<:Space,OT,T} <: $AbstOp{S,OT,T}
            space::S        # the domain space
            order::OT
        end
        immutable $WrappOp{BT<:BandedOperator,S<:Space,OT,T} <: $AbstOp{S,OT,T}
            op::BT
            order::OT
        end


        ## Constructors
        $Op{T}(::Type{T},sp::Space,k)=$Op{typeof(sp),typeof(k),T}(sp,k)
        $Op(::Type{Any},sp::Space,k)=$Op(sp,k)


        $Op(sp::Space{ComplexBasis},k)=$Op{typeof(sp),typeof(k),Complex{real(eltype(domain(sp)))}}(sp,k)
        $Op(sp::Space,k)=$Op{typeof(sp),typeof(k),eltype(domain(sp))}(sp,k)

        $Op(sp::Space)=$Op(sp,1)
        $Op()=$Op(UnsetSpace())
        $Op(k::Number)=$Op(UnsetSpace(),k)

        $Op(d::Domain,n)=$Op(Space(d),n)
        $Op(d::Domain)=$Op(d,1)
        $Op(d::Vector)=$Op(Space(d),1)
        $Op(d::Vector,n)=$Op(Space(d),n)

        function Base.convert{T}(::Type{Operator{T}},D::$Op)
            if T==eltype(D)
                D
            else
                $Op{typeof(D.space),typeof(D.order),T}(D.space,D.order)
            end
        end

        function Base.convert{T}(::Type{BandedOperator{T}},D::$Op)
            if T==eltype(D)
                D
            else
                $Op{typeof(D.space),typeof(D.order),T}(D.space,D.order)
            end
        end

        $WrappOp(op::BandedOperator,order)=$WrappOp{typeof(op),typeof(domainspace(op)),typeof(order),eltype(op)}(op,order)
        $WrappOp(op::BandedOperator)=$WrappOp(op,1)

        function Base.convert{T}(::Type{Operator{T}},D::$WrappOp)
            if T==eltype(D)
                D
            else
                $WrappOp(convert(BandedOperator{T},D.op),D.order)
            end
        end
        function Base.convert{T}(::Type{BandedOperator{T}},D::$WrappOp)
            if T==eltype(D)
                D
            else
                $WrappOp(convert(BandedOperator{T},D.op),D.order)
            end
        end

        ## Routines
        domain(D::$Op)=domain(D.space)
        domainspace(D::$Op)=D.space

        addentries!{OT,T}(::$Op{UnsetSpace,OT,T},A,kr::Range,::Colon)=error("Spaces cannot be inferred for operator")

        function addentries!{S,OT,T}(D::$Op{S,OT,T},A,kr::Range,::Colon)
            # Default is to convert to Canonical and apply operator there
            sp=domainspace(D)
            csp=canonicalspace(sp)
            if conversion_type(csp,sp)==csp   # Conversion(sp,csp) is not banded, or sp==csp
                error("Override addentries! for "*string($Op)*"(::"*string(typeof(sp))*","*string(D.order)*")")
            end
            addentries!(TimesOperator([$Op(csp,D.order),Conversion(sp,csp)]),A,kr,:)
        end

        function bandinds(D::$Op)
            sp=domainspace(D)
            csp=canonicalspace(sp)
            if conversion_type(csp,sp)==csp   # Conversion(sp,csp) is not banded, or sp==csp
                error("Override bandinds for "*string($Op)*"(::"*string(typeof(sp))*","*string(D.order)*")")
            end
            bandinds(TimesOperator([$Op(csp,D.order),Conversion(sp,csp)]))
        end

        # corresponds to default implementation
        function rangespace{S,T}(D::$Op{S,T})
            sp=domainspace(D)
            csp=canonicalspace(sp)
            if conversion_type(csp,sp)==csp   # Conversion(sp,csp) is not banded, or sp==csp
                error("Override rangespace for "*string($Op)*"(::"*string(typeof(sp))*","*string(D.order)*")")
            end
            rangespace($Op(canonicalspace(domainspace(D)),D.order))
        end
        rangespace{T}(D::$Op{UnsetSpace,T})=UnsetSpace()

        #promoting domain space is allowed to change range space
        # for integration, we fall back on existing conversion for now
        promotedomainspace(D::$AbstOp,sp::UnsetSpace)=D
        promotedomainspace(D::$AbstOp,sp::AnySpace)=D


        function promotedomainspace{S<:Space}(D::$AbstOp,sp::S)
            if isambiguous(domain(sp))
                $Op(S(domain(D)),D.order)
            else
                $Op(sp,D.order)
            end
        end

        choosedomainspace(M::$Op{UnsetSpace},sp)=sp  # we assume the space itself will work


        #Wrapper just adds the operator it wraps
        addentries!(D::$WrappOp,A,k::Range,::Colon)=addentries!(D.op,A,k,:)
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







# Multivariate



@calculus_operator(Laplacian)

Laplacian(S::Space,k)=Laplacian{typeof(S),Int,BandedMatrix{eltype(S)}}(S,k)
Laplacian(S)=Laplacian(S,1)
