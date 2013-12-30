

export Domain,tocanonical,fromcanonical,Circle,PeriodicInterval,Interval
export chebyshevpoints


abstract Domain


##General routines

for op = (:tocanonical,:tocanonicalD,:fromcanonical,:fromcanonicalD)
    @eval ($op)(f::AbstractFun,x)=($op)(f.domain,x)
end






## Interval Domains

abstract IntervalDomain <: Domain

chebyshevpoints(n::Integer)= [cos(1.π*k/(n-1)) for k = n-1:-1:0]
points(d::IntervalDomain,n::Integer) =  [fromcanonical(d,cos(1.π*k/(n-1))) for k = n-1:-1:0]
points(d::Vector,n::Integer)=points(Interval(d),n)


####
## Periodic domains

abstract PeriodicDomain <: Domain

points(d::PeriodicDomain,n::Integer) = fromcanonical(d, fourierpoints(n))


fourierpoints(n::Integer)= 1.π*[-1.:2/n:1. - 2/n]


