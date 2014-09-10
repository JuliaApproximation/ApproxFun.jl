

#Ultraspherical Spaces



immutable JacobiWeightSpace{α,β} <: IntervalDomainSpace
    domain::Interval
end

jacobiweight(α,β,x)=(1.+x).^α.*(1.-x).^β
jacobiweight{α,β}(sp::JacobiSpace{α,β})=jacobiweight(α,β,tocanonical(sp,x))

evaluate(sp::JacobiWeightSpace,cfs::Vector,x)=evaluate(ChebyshevSpace(domain(sp)),cfs,x).*jacobiweight(sp,x)
