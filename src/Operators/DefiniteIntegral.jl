export DefiniteIntegral

immutable DefiniteIntegral{D<:FunctionSpace,T<:Number} <: Functional{T}
    domainspace::D
end

# We expect the operator to be real/complex if the basis is real/complex
DefiniteIntegral()=DefiniteIntegral(UnsetSpace())
DefiniteIntegral(dsp::FunctionSpace) = DefiniteIntegral{typeof(dsp),eltype(dsp)}(dsp)
promotedomainspace(::DefiniteIntegral,sp::FunctionSpace)=DefiniteIntegral(sp)

Base.convert{T}(::Type{Functional{T}},Σ::DefiniteIntegral)=DefiniteIntegral{typeof(Σ.domainspace),T}(Σ.domainspace)

domain(Σ::DefiniteIntegral)=domain(Σ.domainspace)
domainspace(Σ::DefiniteIntegral)=Σ.domainspace

getindex(::DefiniteIntegral{UnsetSpace},kr::Range)=error("Spaces cannot be inferred for operator")

