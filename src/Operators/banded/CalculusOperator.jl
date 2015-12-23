export Derivative,Integral,Laplacian


abstract CalculusOperator{S,OT,T}<:BandedOperator{T}


## Note that all functions called in calculus_operator must be exported

iswrapper(::)=false

macro calculus_operator(Op)
    ConcOp=parse("Concrete"*string(Op))
    WrappOp=parse(string(Op)*"Wrapper")
    return esc(quote
        # The SSS, TTT are to work around #9312
        abstract $Op{SSS,OT,TTT} <: CalculusOperator{SSS,OT,TTT}

        immutable $ConcOp{S<:Space,OT,T} <: $Op{S,OT,T}
            space::S        # the domain space
            order::OT
        end
        immutable $WrappOp{BT<:BandedOperator,S<:Space,OT,T} <: $Op{S,OT,T}
            op::BT
            order::OT
        end


        ## Constructors
        $ConcOp(sp::Space,k)=$ConcOp{typeof(sp),typeof(k),promote_type(eltype(sp),eltype(domain(sp)))}(sp,k)
        $Op(sp::Space,k)=$ConcOp(sp,k)

        $Op(sp::Space)=$Op(sp,1)
        $Op()=$Op(UnsetSpace())
        $Op(k::Number)=$Op(UnsetSpace(),k)

        $Op(d::Domain,n)=$Op(Space(d),n)
        $Op(d::Domain)=$Op(d,1)
        $Op(d::Vector)=$Op(Space(d),1)
        $Op(d::Vector,n)=$Op(Space(d),n)
        $ConcOp(S::Space)=$ConcOp(S,1)

        function Base.convert{T}(::Type{Operator{T}},D::$ConcOp)
            if T==eltype(D)
                D
            else
                $ConcOp{typeof(D.space),typeof(D.order),T}(D.space,D.order)
            end
        end

        function Base.convert{T}(::Type{BandedOperator{T}},D::$ConcOp)
            if T==eltype(D)
                D
            else
                $ConcOp{typeof(D.space),typeof(D.order),T}(D.space,D.order)
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
        domain(D::$ConcOp)=domain(D.space)
        domainspace(D::$ConcOp)=D.space

        addentries!{OT,T}(::$ConcOp{UnsetSpace,OT,T},A,kr::Range,::Colon)=error("Spaces cannot be inferred for operator")

        function addentries!{S,OT,T}(D::$ConcOp{S,OT,T},A,kr::Range,::Colon)
            # Default is to convert to Canonical and apply operator there
            sp=domainspace(D)
            csp=canonicalspace(sp)
            if conversion_type(csp,sp)==csp   # Conversion(sp,csp) is not banded, or sp==csp
                error("Override addentries! for "*string($ConcOp)*"(::"*string(typeof(sp))*","*string(D.order)*")")
            end
            addentries!(TimesOperator([$Op(csp,D.order),Conversion(sp,csp)]),A,kr,:)
        end

        function bandinds(D::$ConcOp)
            sp=domainspace(D)
            csp=canonicalspace(sp)
            if conversion_type(csp,sp)==csp   # Conversion(sp,csp) is not banded, or sp==csp
                error("Override bandinds for "*string($ConcOp)*"(::"*string(typeof(sp))*","*string(D.order)*")")
            end
            bandinds(TimesOperator([$Op(csp,D.order),Conversion(sp,csp)]))
        end

        # corresponds to default implementation
        function rangespace{S,T}(D::$ConcOp{S,T})
            sp=domainspace(D)
            csp=canonicalspace(sp)
            if conversion_type(csp,sp)==csp   # Conversion(sp,csp) is not banded, or sp==csp
                error("Override rangespace for "*string($ConcOp)*"(::"*string(typeof(sp))*","*string(D.order)*")")
            end
            rangespace($Op(canonicalspace(domainspace(D)),D.order))
        end
        rangespace{T}(D::$ConcOp{UnsetSpace,T})=UnsetSpace()

        #promoting domain space is allowed to change range space
        # for integration, we fall back on existing conversion for now
        promotedomainspace(D::$Op,sp::UnsetSpace)=D
        promotedomainspace(D::$Op,sp::AnySpace)=D


        function promotedomainspace(D::$Op,sp::Space)
            if isambiguous(domain(sp))
                $Op(typeof(sp)(domain(D)),D.order)
            else
                $Op(sp,D.order)
            end
        end


        #Wrapper just adds the operator it wraps
        addentries!(D::$WrappOp,A,k::Range,::Colon)=addentries!(D.op,A,k,:)
        rangespace(D::$WrappOp)=rangespace(D.op)
        domainspace(D::$WrappOp)=domainspace(D.op)
        bandinds(D::$WrappOp)=bandinds(D.op)

        iswrapper(::$WrappOp)=true
    end)
#     for func in (:rangespace,:domainspace,:bandinds)
#         # We assume the operator wrapped has the correct spaces
#         @eval $func(D::$WrappOp)=$func(D.op)
#     end
end

choosedomainspace(M::CalculusOperator{UnsetSpace},sp)=iswrapper(M)?choosedomainspace(M.op,sp):sp  # we assume the space itself will work



@calculus_operator(Derivative)
@calculus_operator(Integral)

for TYP in (:Derivative,:Integral)
    @eval begin
        function *(D1::$TYP,D2::$TYP)
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
    if typeof(d)==typeof(cd)  || isa(d,PeriodicDomain)
        Integral(space(f))*f
    else
        # map to canonical domain
        setdomain(integrate(setdomain(f,cd)*fromcanonicalD(f)),d)
    end
end

function Base.sum(f::Fun)
    d=domain(f)
    cd=canonicaldomain(d)
    if typeof(cd)==typeof(d)  || isa(d,PeriodicDomain)
        last(cumsum(f))
    else
        # map first
        sum(setdomain(f,cd)*fromcanonicalD(f))
    end
end

function linesum(f::Fun)
    cd=canonicaldomain(f)
    if typeof(cd)==typeof(domain(f))  || isa(d,PeriodicDomain)
        error("override linesum for $(f.space)")
    else
        # map first
        linesum(setdomain(f,cd)*abs(fromcanonicalD(f)))
    end
end




# Multivariate



@calculus_operator(Laplacian)

Laplacian(S::Space,k)=ConcreteLaplacian{typeof(S),Int,BandedMatrix{eltype(S)}}(S,k)
Laplacian(S)=Laplacian(S,1)




## Map to canonical
function defaultderivative(S::Space,order::Integer)
    if typeof(canonicaldomain(S)).name==typeof(domain(S)).name
        # we assume the canonical domain case is implemented
        ConcreteDerivative(S,order)
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
    if typeof(canonicaldomain(sp)).name==typeof(domain(sp)).name
        # we assume the canonical domain case is implemented
        ConcreteIntegral(sp,k)
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
