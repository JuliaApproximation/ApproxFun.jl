export Disk

immutable Disk <: BivariateDomain
    radius::Float64
    center::(Float64,Float64)
end
Disk()=Disk(0.,(0.,0.))

#canonical is rectangle
# we assume radius and centre are zero for now
fromcanonical(D::Disk,x,t)=.5*(1-x)*cos(t),.5*(1-x)*sin(t)
tocanonical(D::Disk,x,y)=1-2sqrt(x^2+y^2),atan2(y,x)


function points(d::Disk,n,m)
    ptsx=gaussjacobi(n,1.,0.)[1]
    ptst=points(PeriodicInterval(),m)
    
    [fromcanonical(d,x,t)[1] for x in ptsx, t in ptst],[fromcanonical(d,x,t)[2] for x in ptsx, t in ptst]
end

immutable DiskSpace <: AbstractProductSpace{JacobiWeightSpace{JacobiSpace},LaurentSpace}
    domain::Disk
end

domain(d::DiskSpace)=d.domain
function space(d::DiskSpace,k::Integer)
    @assert k==2
    LaurentSpace()
end

Space(D::Disk)=DiskSpace(D)

points(d::DiskSpace,n,m)=points(domain(d),n,m)


columnspace(D::DiskSpace,k)=(m=1.div(k,2);JacobiWeightSpace(0.,m,JacobiSpace(2m+1,0.)))

transform(S::DiskSpace,V::Matrix)=transform([columnspace(S,k) for k=1:size(V,2)],LaurentSpace(),V)