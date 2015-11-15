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
        $Op(sp::Space,k)=$Op{typeof(sp),typeof(k),promote_type(eltype(sp),eltype(domain(sp)))}(sp,k)

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
differentiate(f::Fun)=Derivative(space(f))*f
function integrate(f::Fun)
    d=domain(f)
    cd=canonicaldomain(d)
    if typeof(d)==typeof(cd)
        Integral(space(f))*f
    else
        # map to canonical domain
        fc=Fun(f.coefficients,setdomain(space(f),cd))
        x=Fun(identity,cd)
        Mp=fromcanonicalD(f,x)
        g=integrate(fc*Mp)
        Fun(g.coefficients,setdomain(space(g),d))
    end
end

function Base.sum(f::Fun)
    if typeof(canonicaldomain(f))==typeof(domain(f))
        last(cumsum(f))
    else
        # map first
        fc=Fun(f.coefficients,setdomain(space(f),canonicaldomain(f)))
        x=Fun(identity,domain(fc))
        Mp=fromcanonicalD(f,x)
        sum(fc*Mp)
    end
end

function linesum(f::Fun)
    if typeof(canonicaldomain(f))==typeof(domain(f))
        error("override linesum for $(f.space)")
    else
        # map first
        fc=Fun(f.coefficients,setdomain(space(f),canonicaldomain(f)))
        x=Fun(identity,domain(fc))
        Mp=fromcanonicalD(f,x)
        linesum(fc*abs(Mp))
    end
end




# Multivariate



@calculus_operator(Laplacian)

Laplacian(S::Space,k)=Laplacian{typeof(S),Int,BandedMatrix{eltype(S)}}(S,k)
Laplacian(S)=Laplacian(S,1)




## Map to canonical
function defaultderivative(S::Space,order::Integer)
    if typeof(canonicaldomain(S))==typeof(domain(S))
        # we assume the canonical domain case is implemented
        Derivative{typeof(S),typeof(order),promote_type(eltype(S),eltype(domain(S)))}(S,order)
    else
        D1=invfromcanonicalD(S)*Derivative(setdomain(S,canonicaldomain(S)))
        D=DerivativeWrapper(SpaceOperator(D1,S,setdomain(rangespace(D1),domain(S))),1)
        if order==1
            D
        else
            DerivativeWrapper(TimesOperator(Derivative(rangespace(D),order-1),D),order)
        end
    end
end



Derivative(S::Space,order::Integer)=defaultderivative(S,order)


function Integral(sp::Space,k::Integer)
    if typeof(canonicaldomain(sp))==typeof(domain(sp))
        # we assume the canonical domain case is implemented
        Integral{typeof(sp),typeof(k),promote_type(eltype(sp),eltype(domain(sp)))}(sp,k)
    elseif k > 1
        Q=Integral(sp,1)
        IntegralWrapper(TimesOperator(Integral(rangespace(Q),k-1),Q),k)
    else # k==1
        csp=setdomain(sp,canonicaldomain(sp))

        x=Fun(identity,domain(csp))
        M=Multiplication(fromcanonicalD(sp,x),csp)
        Q=Integral(rangespace(M))*M
        IntegralWrapper(SpaceOperator(Q,sp,setdomain(rangespace(Q),domain(sp))),1)
    end
end
