##
# LogWeight
# represents a function on [-1,1] weighted by log((1+x)^α*(1-x)^β)
##


immutable LogWeight{S<:IntervalSpace} <: WeightSpace
    α::Float64
    β::Float64
    space::S
end

spacescompatible(A::LogWeight,B::LogWeight)=A.α==B.α && A.β == B.β && spacescompatible(A.space,B.space)

logweight(α,β,x)=log((1+x).^α.*(1-x).^β)
weight(sp::LogWeight,x)=logweight(sp.α,sp.β,tocanonical(sp,x))


setdomain(sp::LogWeight,d::Domain)=LogWeight(sp.α,sp.β,setdomain(sp.space,d))

function coefficients(f::Vector,sp1::LogWeight,sp2::LogWeight)
    α,β=sp1.α,sp1.β
    c,d=sp2.α,sp2.β

    if isapprox(c,α) && isapprox(d,β)
        coefficients(f,sp1.space,sp2.space)
    else
        (Conversion(sp1,sp2)*f)
    end
end

