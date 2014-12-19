typealias PeriodicSpace{T} FunctionSpace{T,PeriodicInterval}       
canonicaldomain{T<:PeriodicSpace}(::Type{T})=PeriodicInterval()

## Toeplitz

ToeplitzOperator{T,D<:PeriodicSpace}(f::Fun{D,T})=ToeplitzOperator(f.coefficients|>deinterlace)
LaurentOperator{T,D<:PeriodicSpace}(f::Fun{D,T})=LaurentOperator(flipud(f.coefficients|>deinterlace))

