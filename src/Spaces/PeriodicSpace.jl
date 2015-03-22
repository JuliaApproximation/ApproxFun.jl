typealias PeriodicSpace{T} FunctionSpace{T,PeriodicInterval}
canonicaldomain{T<:PeriodicSpace}(::Type{T})=PeriodicInterval()

## Evaluation

Evaluation(d::PeriodicDomain,x::Number,n...)=Evaluation(Laurent(d),complex(x),n...)

## Definite Integral

DefiniteIntegral(d::PeriodicDomain)=DefiniteIntegral(Laurent(d))
DefiniteLineIntegral(d::PeriodicDomain)=DefiniteLineIntegral(Laurent(d))

## Toeplitz

ToeplitzOperator{T,D<:PeriodicSpace}(f::Fun{D,T})=ToeplitzOperator(f.coefficients[2:2:end],f.coefficients[1:2:end])
LaurentOperator{T,D<:PeriodicSpace}(f::Fun{D,T})=LaurentOperator(f.coefficients[3:2:end],f.coefficients[[1;2:2:end]])

