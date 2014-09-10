

#Ultraspherical Spaces



immutable JacobiWeightSpace{α,β} <: IntervalDomainSpace
    domain::Interval
end

jacobiweight(α,β,x)=(1.+x).^α.*(1.-x).^β
