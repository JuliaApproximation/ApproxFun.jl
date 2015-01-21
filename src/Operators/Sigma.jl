export Σ
# The Σ operator is intended to act as the definite integration
# operator, expressing the operation of integrating an unknown Fun with
# a bivariate kernel. For orthogonal cases, the operator is banded.

immutable Σ{T<:Number,D<:FunctionSpace,R<:FunctionSpace} <: BandedOperator{T}
    domainspace::D
    rangespace::R
end

# We expect the operator to be real/complex if the basis is real/complex
Σ()=Σ(AnySpace(),AnySpace())
Σ(dsp::FunctionSpace,rsp::FunctionSpace) = Σ{promote_type(eltype(dsp),eltype(rsp)),typeof(dsp),typeof(rsp)}(dsp,rsp)

Base.convert{T}(::Type{BandedOperator{T}},S::Σ)=Σ{T,typeof(S.domainspace),typeof(S.rangespace)}(S.domainspace,S.rangespace)

domain(S::Σ)=domain(S.domainspace)
domainspace(S::Σ)=S.domainspace
rangespace(S::Σ)=S.rangespace

addentries!{T}(::Σ{T,AnySpace,AnySpace},A,kr::Range)=error("Spaces cannot be inferred for operator")

bandinds(S::Σ) = 0,0
