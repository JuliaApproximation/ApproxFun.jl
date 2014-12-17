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

Σ(d::PeriodicDomain)=Σ(LaurentSpace(d),LaurentSpace(d))
Σ(d::IntervalDomain)=Σ(JacobiWeightSpace(-.5,-.5,ChebyshevSpace(d)),ChebyshevSpace(d))
function Σ(α::Number,β::Number,d::IntervalDomain)
    @assert α == β
    @assert int(α+.5) == α+.5
    Σ(JacobiWeightSpace(α,β,UltrasphericalSpace{int(α+.5)}(d)),UltrasphericalSpace{int(α+.5)}(d))
end
Σ(α::Number,β::Number) = Σ(α,β,Interval())

domain(S::Σ)=domain(S.domainspace)
domainspace(S::Σ)=S.domainspace
rangespace(S::Σ)=S.rangespace

addentries!{T}(::Σ{T,AnySpace,AnySpace},A::ShiftArray,kr::Range)=error("Spaces cannot be inferred for operator")

bandinds(S::Σ) = 0,0