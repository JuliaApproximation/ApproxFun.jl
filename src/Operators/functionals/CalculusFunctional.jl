export DefiniteIntegral,DefiniteLineIntegral

abstract type CalculusFunctional{S,T} <: Operator{T} end

@functional CalculusFunctional

##TODO: Add ConcreteOp

macro calculus_functional(Op)
    ConcOp=Meta.parse("Concrete"*string(Op))
    WrappOp=Meta.parse(string(Op)*"Wrapper")
    return esc(quote
        abstract type $Op{SSS,TTT} <: CalculusFunctional{SSS,TTT} end
        struct $ConcOp{S,T} <: $Op{S,T}
            domainspace::S
        end
        struct $WrappOp{BT<:Operator,S<:Space,T} <: $Op{S,T}
            op::BT
        end

        @wrapper $WrappOp


        # We expect the operator to be real/complex if the basis is real/complex
        $ConcOp(dsp::Space) = $ConcOp{typeof(dsp),prectype(dsp)}(dsp)

        $Op() = $Op(UnsetSpace())
        $Op(dsp) = $ConcOp(dsp)
        $Op(d::Domain) = $Op(Space(d))

        ApproxFun.promotedomainspace(::$Op,sp::Space) = $Op(sp)


        Base.convert(::Type{Operator{T}},Σ::$ConcOp) where {T} =
            (T==eltype(Σ) ? Σ : $ConcOp{typeof(Σ.domainspace),T}(Σ.domainspace))::Operator{T}

        ApproxFun.domain(Σ::$ConcOp) = domain(Σ.domainspace)
        ApproxFun.domainspace(Σ::$ConcOp) = Σ.domainspace

        Base.getindex(::$ConcOp{UnsetSpace},kr::AbstractRange) =
            error("Spaces cannot be inferred for operator")

        $WrappOp(op::Operator) =
            $WrappOp{typeof(op),typeof(domainspace(op)),eltype(op)}(op)


        Base.convert(::Type{Operator{T}},Σ::$WrappOp) where {T} =
            (T==eltype(Σ) ? Σ : $WrappOp(convert(Operator{T},Σ.op)))::Operator{T}
    end)
end

@calculus_functional(DefiniteIntegral)
@calculus_functional(DefiniteLineIntegral)


#default implementation

DefiniteIntegral(sp::UnsetSpace) = ConcreteDefiniteIntegral(sp)
DefiniteLineIntegral(sp::UnsetSpace) = ConcreteDefiniteLineIntegral(sp)

function DefiniteIntegral(sp::Space)
    if typeof(canonicaldomain(sp)) == typeof(domain(sp))
        # try using `Integral`
        Q = Integral(sp)
        rsp = rangespace(Q)
        DefiniteIntegralWrapper((Evaluation(rsp,last)-Evaluation(rsp,first))*Q)
    else
        # try mapping to canonical domain
        M = Multiplication(fromcanonicalD(sp),setcanonicaldomain(sp))
        Op = DefiniteIntegral(rangespace(M))*M
        DefiniteIntegralWrapper(SpaceOperator(Op,sp,rangespace(Op)))
    end
end

function DefiniteLineIntegral(sp::Space)
    if typeof(canonicaldomain(sp)) == typeof(domain(sp))
        error("Override DefiniteLineIntegral for $sp")
    end

    M = Multiplication(abs(fromcanonicalD(sp)),setcanonicaldomain(sp))
    Op = DefiniteLineIntegral(rangespace(M))*M
    DefiniteLineIntegralWrapper(SpaceOperator(Op,sp,rangespace(Op)))
end


#TODO: Remove SPECIALOPS reimplement
# *{T,D<:DefiniteIntegral,M<:Multiplication}(A::TimesFunctional{T,D,M},b::Fun) = bilinearform(A.op.f,b)
# *{T,D<:DefiniteLineIntegral,M<:Multiplication}(A::TimesFunctional{T,D,M},b::Fun) = linebilinearform(A.op.f,b)
# *{T,D<:Union{DefiniteIntegral,DefiniteLineIntegral},
#   M<:Multiplication,V}(A::FunctionalOperator{TimesFunctional{T,D,M},V},b::Fun) =
#     Fun(A.op*b)
