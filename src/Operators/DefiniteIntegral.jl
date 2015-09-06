export DefiniteIntegral,DefiniteLineIntegral

abstract CalculusFunctional{S,T} <: Functional{T}

macro calculus_functional(Func)
    AbstOp=parse("Abstract"*string(Func))
    WrappOp=parse(string(Func)*"Wrapper")
    return esc(quote
        abstract $AbstOp{SSS,TTT} <: CalculusFunctional{SSS,TTT}
        immutable $Func{S,T} <: $AbstOp{S,T}
            domainspace::S
        end
        immutable $WrappOp{BT<:Functional,S<:FunctionSpace,T} <: $AbstOp{S,T}
            func::BT
        end


        # We expect the operator to be real/complex if the basis is real/complex
        $Func()=$Func(UnsetSpace())
        $Func(dsp::FunctionSpace) = $Func{typeof(dsp),eltype(dsp)}(dsp)
        $Func(d::Domain) = $Func(Space(d))

        promotedomainspace(::$AbstOp,sp::FunctionSpace)=$Func(sp)

        Base.convert{OT<:$Func}(::Type{OT},Σ::OT)=Σ
        Base.convert{OT<:Operator}(::Type{OT},Σ::$Func)=$Func{typeof(Σ.domainspace),eltype(OT)}(Σ.domainspace)::OT

        domain(Σ::$Func)=domain(Σ.domainspace)
        domainspace(Σ::$Func)=Σ.domainspace

        getindex(::$Func{UnsetSpace},kr::Range)=error("Spaces cannot be inferred for operator")

        $WrappOp(op::Functional)=$WrappOp{typeof(op),typeof(domainspace(op)),eltype(op)}(op)

        Base.convert{OT<:$WrappOp}(::Type{OT},Σ::OT)=Σ
        Base.convert{OT<:Operator}(::Type{OT},Σ::$WrappOp)=$WrappOp(convert(Operator{eltype(OT)},Σ.func))::OT

        #Wrapper just adds the operator it wraps
        getindex(D::$WrappOp,k::Range)=D.func[k]
        domainspace(D::$WrappOp)=domainspace(D.func)
        datalength(D::$WrappOp)=datalength(D.func)
    end)
end

@calculus_functional(DefiniteIntegral)
@calculus_functional(DefiniteLineIntegral)


#default implementation
function getindex(B::DefiniteIntegral,kr::Range)
    S=domainspace(B)
    Q=Integral(S)
    A=(Evaluation(S,true)-Evaluation(S,false))*Q
    A[kr]
end
