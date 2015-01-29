export Σ
# The Σ operator is intended to act as the definite integration
# operator, expressing the operation of integrating an unknown Fun with
# a bivariate kernel. For orthogonal cases, the operator is banded.

immutable Σ{D<:FunctionSpace,T<:Number} <: Functional{T}
    domainspace::D
end

typealias DefiniteIntegral Σ

# We expect the operator to be real/complex if the basis is real/complex
Σ()=Σ(UnsetSpace())
Σ(dsp::FunctionSpace) = Σ{typeof(dsp),promote_type(eltype(dsp))}(dsp)
promotedomainspace(::Σ,sp::FunctionSpace)=Σ(sp)

Base.convert{T}(::Type{Functional{T}},S::Σ)=Σ{typeof(S.domainspace),T}(S.domainspace)

domain(S::Σ)=domain(S.domainspace)
domainspace(S::Σ)=S.domainspace

getindex(::Σ{UnsetSpace},kr::Range)=error("Spaces cannot be inferred for operator")

