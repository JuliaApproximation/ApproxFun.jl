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
Σ{T1,T2}(dsp::FunctionSpace{T1},rsp::FunctionSpace{T2}) = Σ{promote_type(T1,T2),typeof(dsp),typeof(rsp)}(dsp,rsp)



domain(S::Σ)=domain(S.domainspace)
domainspace(S::Σ)=S.domainspace
rangespace(S::Σ)=S.rangespace

addentries!{T}(::Σ{T,AnySpace,AnySpace},A,kr::Range)=error("Spaces cannot be inferred for operator")

bandinds(S::Σ) = 0,0