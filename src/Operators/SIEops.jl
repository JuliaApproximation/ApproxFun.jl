export Σ

immutable Σ{S<:FunctionSpace,T<:Number} <: BandedOperator{T}
    space::S        # the domain space
end

Σ{S<:PeriodicDomainSpace}(sp::S) = Σ{S,Complex{Float64}}(sp)
Σ{S<:FunctionSpace}(sp::S)=Σ{S,Float64}(sp)
Σ()=Σ(AnySpace())

Σ(d::PeriodicDomain)=Σ(LaurentSpace(d))
Σ(d::IntervalDomain)=Σ(ChebyshevSpace(d))

domain(S::Σ)=domain(S.space)
domainspace(S::Σ)=S.space
rangespace(S::Σ)=S.space

addentries!{T}(::Σ{AnySpace,T},A::ShiftArray,kr::Range)=error("Spaces cannot be inferred for operator")

bandinds(S::Σ) = 0,0