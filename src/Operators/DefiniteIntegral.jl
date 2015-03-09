export DefiniteIntegral,DefiniteLineIntegral

macro calculus_functional(Func)
    return esc(quote
        immutable $Func{S<:FunctionSpace,T<:Number} <: Functional{T}
            domainspace::S
        end

        # We expect the operator to be real/complex if the basis is real/complex
        $Func()=$Func(UnsetSpace())
        $Func(dsp::FunctionSpace) = $Func{typeof(dsp),eltype(dsp)}(dsp)

        promotedomainspace(::$Func,sp::FunctionSpace)=$Func(sp)

        Base.convert{T}(::Type{Functional{T}},Σ::$Func)=$Func{typeof(Σ.domainspace),T}(Σ.domainspace)

        domain(Σ::$Func)=domain(Σ.domainspace)
        domainspace(Σ::$Func)=Σ.domainspace

        getindex(::$Func{UnsetSpace},kr::Range)=error("Spaces cannot be inferred for operator")

    end)
end

@calculus_functional(DefiniteIntegral)
@calculus_functional(DefiniteLineIntegral)
