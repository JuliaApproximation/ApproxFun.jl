

immutable Disk <: BivariateDomain
    radius::Float64
    center::(Float64,Float64)
end
Disk()=Disk(0.,(0.,0.))

#canonical is rectangle
# we assume radius and centre are zero for now
fromcanonical(D::Disk,x,t)=.5*(1-x)*cos(t),.5*(1-x)*sin(t)
tocanonical(D::Disk,x,y)=1-2sqrt(x^2+y^2),atan2(y,x)


immutable DiskSpace <: AbstractProductSpace{JacobiSpace,LaurentSpace}
    domain::Disk
end

domain(d::DiskSpace)=d.domain


