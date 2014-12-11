abstract PeriodicDomainSpace{T} <: DomainSpace{T,PeriodicInterval}       


## Toeplitz

ToeplitzOperator{T,D<:PeriodicDomainSpace}(f::Fun{D,T})=ToeplitzOperator(f.coefficients|>deinterlace)
LaurentOperator{T,D<:PeriodicDomainSpace}(f::Fun{D,T})=LaurentOperator(flipud(f.coefficients|>deinterlace))

