export Σ
# The Σ operator is intended to act as the definite integration
# operator, expressing the operation of integrating an unknown Fun with
# a bivariate kernel. For orthogonal cases, the operator is banded.

immutable Σ{T<:Number,D<:FunctionSpace} <: Functional{T}
    domainspace::D
end

# We expect the operator to be real/complex if the basis is real/complex
Σ()=Σ(UnsetSpace())
Σ(dsp::FunctionSpace) = Σ{promote_type(eltype(dsp)),typeof(dsp)}(dsp)

Base.convert{T}(::Type{Functional{T}},S::Σ)=Σ{T,typeof(S.domainspace)}(S.domainspace)

domain(S::Σ)=domain(S.domainspace)
domainspace(S::Σ)=S.domainspace

getindex{T}(::Σ{T,UnsetSpace},kr::Range)=error("Spaces cannot be inferred for operator")

