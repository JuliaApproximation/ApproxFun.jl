

immutable FourierSpace{T<:Union(PeriodicDomain,DataType)} <: PeriodicDomainSpace
    domain::T
end

FourierSpace()=FourierSpace(Any)

==(a::FourierSpace,b::FourierSpace)= a.domain==b.domain

