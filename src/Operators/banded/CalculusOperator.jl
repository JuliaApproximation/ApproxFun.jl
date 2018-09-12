export Derivative,Integral,Laplacian,Volterra


abstract type CalculusOperator{S,OT,T}<:Operator{T} end


## Note that all functions called in calculus_operator must be exported

macro calculus_operator(Op)
    ConcOp=Meta.parse("Concrete"*string(Op))
    WrappOp=Meta.parse(string(Op)*"Wrapper")
    DefaultOp=Meta.parse("Default"*string(Op))
    return esc(quote
        # The SSS, TTT are to work around #9312
        abstract type $Op{SSS,OT,TTT} <: CalculusOperator{SSS,OT,TTT} end

        struct $ConcOp{S<:Space,OT,T} <: $Op{S,OT,T}
            space::S        # the domain space
            order::OT
        end
        struct $WrappOp{BT<:Operator,S<:Space,OT,T} <: $Op{S,OT,T}
            op::BT
            order::OT
        end

        @wrapper $WrappOp


        ## Constructors
        $ConcOp(sp::Space,k) = $ConcOp{typeof(sp),typeof(k),prectype(sp)}(sp,k)

        $Op(sp::UnsetSpace,k) = $ConcOp(sp,k)
        $Op(sp::UnsetSpace,k::Real) = $ConcOp(sp,k)
        $Op(sp::UnsetSpace,k::Integer) = $ConcOp(sp,k)

        function $DefaultOp(sp::Space,k)
            csp=canonicalspace(sp)
            if conversion_type(csp,sp)==csp   # Conversion(sp,csp) is not banded, or sp==csp
               error("Implement $(string($Op))($(string(sp)),$k)")
            end
            $WrappOp(TimesOperator([$Op(csp,k),Conversion(sp,csp)]),k)
        end

        $DefaultOp(d,k) = $Op(Space(d),k)

        $DefaultOp(sp) = $Op(sp,1)
        $DefaultOp() = $Op(UnsetSpace())
        $DefaultOp(k::Number) = $Op(UnsetSpace(),k)
        $DefaultOp(k::Vector) = $Op(UnsetSpace(),k)

        $Op(x...) = $DefaultOp(x...)
        $ConcOp(S::Space) = $ConcOp(S,1)

        function Base.convert(::Type{Operator{T}},D::$ConcOp) where T
            if T==eltype(D)
                D
            else
                $ConcOp{typeof(D.space),typeof(D.order),T}(D.space,D.order)
            end
        end

        $WrappOp(op::Operator,order) =
            $WrappOp{typeof(op),typeof(domainspace(op)),typeof(order),eltype(op)}(op,order)
        $WrappOp(op::Operator) = $WrappOp(op,1)

        function Base.convert(::Type{Operator{T}},D::$WrappOp) where T
            if T==eltype(D)
                D
            else
                # work around typeinfernece bug
                op=convert(Operator{T},D.op)
                $WrappOp{typeof(op),typeof(domainspace(op)),typeof(D.order),T}(op,D.order)
            end
        end

        ## Routines
        ApproxFun.domain(D::$ConcOp) = domain(D.space)
        ApproxFun.domainspace(D::$ConcOp) = D.space

        Base.getindex(::$ConcOp{UnsetSpace,OT,T},k::Integer,j::Integer) where {OT,T} =
            error("Spaces cannot be inferred for operator")
        ApproxFun.rangespace(D::$ConcOp{UnsetSpace,T}) where {T} = UnsetSpace()

        #promoting domain space is allowed to change range space
        # for integration, we fall back on existing conversion for now
        ApproxFun.promotedomainspace(D::$Op,sp::UnsetSpace) = D


        function ApproxFun.promotedomainspace(D::$Op,sp::Space)
            if isambiguous(domain(sp))
                $Op(typeof(sp)(domain(D)),D.order)
            else
                $Op(sp,D.order)
            end
        end
    end)
#     for func in (:rangespace,:domainspace,:bandwidths)
#         # We assume the operator wrapped has the correct spaces
#         @eval $func(D::$WrappOp)=$func(D.op)
#     end
end

choosedomainspace(M::CalculusOperator{UnsetSpace},sp::Space) =
    iswrapper(M) ? choosedomainspace(M.op,sp) : sp  # we assume the space itself will work



@calculus_operator(Derivative)
@calculus_operator(Integral)
@calculus_operator(Volterra)








## Overrideable


## Convenience routines


integrate(d::Domain) = Integral(d,1)


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
        g=integrate(f)
        last(g)-first(g)
    else
        # map first
        sum(setdomain(f,cd)*fromcanonicalD(f))
    end
end

function linesum(f::Fun)
    cd=canonicaldomain(f)
    d=domain(f)

    if isreal(d)
        a,b=first(d),last(d)
        sign(last(b)-first(a))*sum(f)
    elseif typeof(cd)==typeof(d)  || isa(d,PeriodicDomain)
        error("override linesum for $(f.space)")
    else
        # map first
        linesum(setdomain(f,cd)*abs(fromcanonicalD(f)))
    end
end




# Multivariate



@calculus_operator(Laplacian)


## Map to canonical
function DefaultDerivative(sp::Space,k::Integer)
    if typeof(canonicaldomain(sp)).name==typeof(domain(sp)).name
        # this is the normal default constructor
        csp=canonicalspace(sp)
        if conversion_type(csp,sp)==csp   # Conversion(sp,csp) is not banded, or sp==csp
           error("Implement Derivative($(string(sp)),$k)")
        end
        DerivativeWrapper(TimesOperator([Derivative(csp,k),Conversion(sp,csp)]),k)
    else
        D1=invfromcanonicalD(sp)*Derivative(setdomain(sp,canonicaldomain(sp)))
        D=DerivativeWrapper(SpaceOperator(D1,sp,setdomain(rangespace(D1),domain(sp))),1)
        if k==1
            D
        else
            DerivativeWrapper(TimesOperator(Derivative(rangespace(D),k-1),D),k)
        end
    end
end


function DefaultIntegral(sp::Space,k::Integer)
    if typeof(canonicaldomain(sp)).name==typeof(domain(sp)).name
        # this is the normal default constructor
        csp=canonicalspace(sp)
        if conversion_type(csp,sp)==csp   # Conversion(sp,csp) is not banded, or sp==csp
            # we require that Integral is overridden
            error("Implement Integral($(string(sp)),$k)")
        end
        IntegralWrapper(TimesOperator([Integral(csp,k),Conversion(sp,csp)]),k)
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


for TYP in (:Derivative,:Integral,:Laplacian)
    @eval begin
        function *(D1::$TYP,D2::$TYP)
            @assert domain(D1) == domain(D2)

            $TYP(domainspace(D2),D1.order+D2.order)
        end
    end
end


"""
`Derivative(sp::Space,k::Int)` represents the `k`-th derivative on `sp`.
"""
Derivative(::Space,::Int)

"""
`Derivative(sp::Space,k::Vector{Int})` represents a partial derivative on a multivariate space.
For example,
```julia
Dx = Derivative(Chebyshev()^2,[1,0]) # ∂/∂x
Dy = Derivative(Chebyshev()^2,[0,1]) # ∂/∂y
```
"""
Derivative(::Space,::Vector{Int})

"""
`Derivative(sp::Space)` represents the first derivative `Derivative(sp,1)`.
"""
Derivative(::Space)

"""
`Derivative(k)` represents the `k`-th derivative, acting on an unset space.
Spaces will be inferred when applying or manipulating the operator.
"""
Derivative(k)

"""
`Derivative()` represents the first derivative on an unset space.
Spaces will be inferred when applying or manipulating the operator.
"""
Derivative()


"""
`Integral(sp::Space,k::Int)` represents a `k`-th integral on `sp`.
There is no guarantee on normalization.
"""
Integral(::Space,::Int)


"""
`Integral(sp::Space)` represents the first integral `Integral(sp,1)`.
"""
Integral(::Space)

"""
Integral(k)` represents the `k`-th integral, acting on an unset space.
Spaces will be inferred when applying or manipulating the operator.
"""
Integral(k)

"""
`Intergral()` represents the first integral on an unset space.
Spaces will be inferred when applying or manipulating the operator.
"""
Integral()

"""
`Laplacian(sp::Space)` represents the laplacian on space `sp`.
"""
Laplacian(::Space)

"""
`Laplacian()` represents the laplacian on an unset space.
Spaces will be inferred when applying or manipulating the operator.
"""
Laplacian()
