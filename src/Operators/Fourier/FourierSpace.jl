

type FourierSpace{T<:Union(PeriodicDomain,DataType)} <: OperatorSpace
    domain::T
end

FourierSpace()=FourierSpace(Any)

==(a::FourierSpace,b::FourierSpace)= a.domain==b.domain

