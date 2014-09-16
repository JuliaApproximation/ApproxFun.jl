export EvaluationFunctional

## EvaluationFunctional constructors

type EvaluationFunctional{T<:Number,S<:FunctionSpace} <: Functional{T}
    space::S
    x::T
    order::Int
end




EvaluationFunctional(d::FunctionSpace,x::Number)=EvaluationFunctional(d,x,0)
EvaluationFunctional(d::IntervalDomain,x::Number,n...)=EvaluationFunctional(ChebyshevSpace(d),x,n...)
EvaluationFunctional(d::PeriodicDomain,x::Number,n...)=EvaluationFunctional(LaurentSpace(d),complex(x),n...)
EvaluationFunctional{T<:Number}(d::Vector{T},x::Number,o::Integer)=EvaluationFunctional(Interval(d),x,o)
EvaluationFunctional(x::Number)=EvaluationFunctional(Interval(),x,0)

domainspace(E::EvaluationFunctional)=E.space
domain(E::EvaluationFunctional)=domain(E.space)