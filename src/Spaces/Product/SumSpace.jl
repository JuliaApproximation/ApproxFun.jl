

immutable PeriodicSumSpace{S<:DomainSpace,T<:DomainSpace} <: PeriodicDomainSpace
    spaces::(S,T)
end

function PeriodicSumSpace{S<:DomainSpace,T<:DomainSpace}(A::S,B::T)
    @assert domain(A)==domain(B)
    PeriodicSumSpace{S,T}((A,B))
end

##TODO IntervalSumSpace
typealias SumSpace PeriodicSumSpace




domain(A::SumSpace)=domain(A.spaces[1])
evaluate{T,D<:SumSpace}(f::Fun{T,D},x)=evaluate(Fun(f.coefficients[1:2:end],space(f).spaces[1]),x)+evaluate(Fun(f.coefficients[2:2:end],space(f).spaces[2]),x)


=={S,T}(A::SumSpace{S,T},B::SumSpace{S,T})=A.spaces[1]==B.spaces[1] && A.spaces[2]==B.spaces[2]
