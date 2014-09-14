

immutable SumSpace{S<:DomainSpace,T<:DomainSpace} <: DomainSpace
    spaces::(S,T)
end

function SumSpace{S<:DomainSpace,T<:DomainSpace}(A::S,B::T)
    @assert domain(A)==domain(B)
    SumSpace{S,T}((A,B))
end


domain(A::SumSpace)=domain(A.spaces[1])
evaluate{T,D<:SumSpace}(f::IFun{T,D},x)=evaluate(IFun(f.coefficients[1:2:end],space(f).spaces[1]),x)+evaluate(IFun(f.coefficients[2:2:end],space(f).spaces[2]),x)


=={S,T}(A::SumSpace{S,T},B::SumSpace{S,T})=A.spaces[1]==B.spaces[1] && A.spaces[2]==B.spaces[2]
