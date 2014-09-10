

#Ultraspherical Spaces



immutable JacobiWeightSpace{α,β} <: IntervalDomainSpace
    domain::Interval
end

jacobiweight(α,β,x)=(1.+x).^α.*(1.-x).^β
jacobiweight{α,β}(sp::JacobiSpace{α,β},x)=jacobiweight(α,β,tocanonical(sp,x))

evaluate{T,J<:JacobiWeightSpace}(f::IFun{T,J},x)=jacobiweight(space(f),x).*IFun(f.coefficients,domain(f))[x]
