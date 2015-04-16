abstract PeriodicSpace{T} <: FunctionSpace{T}
canonicaldomain{T<:PeriodicSpace}(::Type{T})=PeriodicInterval()

## Evaluation

Evaluation(d::PeriodicDomain,x::Number,n...)=Evaluation(Laurent(d),complex(x),n...)

## Definite Integral

DefiniteIntegral(d::PeriodicDomain)=DefiniteIntegral(Laurent(d))
DefiniteLineIntegral(d::PeriodicDomain)=DefiniteLineIntegral(Laurent(d))

## Toeplitz

