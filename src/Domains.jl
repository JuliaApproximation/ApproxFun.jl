

export Domain,tocanonical,fromcanonical,Circle,PeriodicInterval,Interval
export chebyshevpoints


abstract Domain


##General routines

for op = (:tocanonical,:tocanonicalD,:fromcanonical,:fromcanonicalD)
    @eval ($op)(f::AbstractFun,x)=($op)(f.domain,x)
end






## Interval Domains

abstract IntervalDomain <: Domain

chebyshevpoints(n::Integer)= cos(1.π*[n-1:-1:0]/(n-1))
points(d::IntervalDomain,n::Integer) = fromcanonical(d, chebyshevpoints(n))



####
## Periodic domains

abstract PeriodicDomain <: Domain

points(d::PeriodicDomain,n::Integer) = fromcanonical(d, fourierpoints(n))


fourierpoints(n::Integer)= 1.π*[-1.:2/n:1. - 2/n]


