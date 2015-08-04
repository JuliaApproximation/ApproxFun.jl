export DefiniteIntegral,DefiniteLineIntegral

macro calculus_functional(Func,AbstFunc,WrappFun)
    return esc(quote
        immutable $Func{S,T} <: Functional{T}
            domainspace::S
        end

        # We expect the operator to be real/complex if the basis is real/complex
        $Func()=$Func(UnsetSpace())
        $Func(dsp::FunctionSpace) = $Func{typeof(dsp),eltype(dsp)}(dsp)
        $Func(d::Domain) = $Func(Space(d))

        promotedomainspace(::$Func,sp::FunctionSpace)=$Func(sp)

        Base.convert{OT<:$Func}(::Type{OT},Σ::OT)=Σ
        Base.convert{OT<:Operator}(::Type{OT},Σ::$Func)=$Func{typeof(Σ.domainspace),eltype(OT)}(Σ.domainspace)::OT

        domain(Σ::$Func)=domain(Σ.domainspace)
        domainspace(Σ::$Func)=Σ.domainspace

        getindex(::$Func{UnsetSpace},kr::Range)=error("Spaces cannot be inferred for operator")

    end)
end

@calculus_functional(DefiniteIntegral,AbstractDefiniteIntegral,DefiniteIntegralWrapper)
@calculus_functional(DefiniteLineIntegral,AbstractDefiniteLineIntegral,DefiniteLineIntegralWrapper)


#default implementation
function getindex(B::DefiniteIntegral,kr::Range)
    S=domainspace(B)
    Q=Integral(S)
    A=(Evaluation(S,true)-Evaluation(S,false))*Q
    A[kr]
end
