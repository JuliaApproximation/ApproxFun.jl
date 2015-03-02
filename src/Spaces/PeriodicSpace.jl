typealias PeriodicSpace{T} FunctionSpace{T,PeriodicInterval}
canonicaldomain{T<:PeriodicSpace}(::Type{T})=PeriodicInterval()

## Evaluation

Evaluation(d::PeriodicDomain,x::Number,n...)=Evaluation(Laurent(d),complex(x),n...)

## Definite Integral

DefiniteIntegral(d::PeriodicDomain)=DefiniteIntegral(Laurent(d))

## Toeplitz

ToeplitzOperator{T,D<:PeriodicSpace}(f::Fun{D,T})=ToeplitzOperator(f.coefficients|>deinterlace)
LaurentOperator{T,D<:PeriodicSpace}(f::Fun{D,T})=LaurentOperator(flipud(f.coefficients|>deinterlace))

