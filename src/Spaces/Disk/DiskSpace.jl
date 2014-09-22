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


function points(d::Disk,n,m,k)
    ptsx=gaussjacobi(n,1.,0.)[1]
    ptst=points(PeriodicInterval(),m)
    
    Float64[fromcanonical(d,x,t)[k] for x in ptsx, t in ptst]
end
points(d::BivariateDomain,n,m)=points(d,n,m,1),points(d,n,m,2)

immutable DiskSpace{S<:PeriodicDomainSpace} <: AbstractProductSpace{JacobiWeightSpace{JacobiSpace},LaurentSpace}
    domain::Disk
    spacet::S
end

DiskSpace(D::Disk)=DiskSpace(D,FourierSpace())

domain(d::DiskSpace)=d.domain
function space(D::DiskSpace,k::Integer)
    @assert k==2
    D.spacet
end

Space(D::Disk)=DiskSpace(D)

points(d::DiskSpace,n...)=points(domain(d),n...)


columnspace(D::DiskSpace,k)=(m=1.div(k,2);JacobiWeightSpace(0.,m,JacobiSpace(2m+1,0.)))

transform(S::DiskSpace,V::Matrix)=transform([columnspace(S,k) for k=1:size(V,2)],S.spacet,V)


#Base.real{S,LaurentSpace,T}(u::ProductFun{S,V,T})=real(TensorFun(real(u.coefficients),space(u,2)).').'-imag(TensorFun(imag(u.coefficients),space(u,2)).').'
#Base.imag{S,LaurentSpace,T}(u::ProductFun{S,V,T})=real(TensorFun(imag(u.coefficients),space(u,2)).').'+imag(TensorFun(real(u.coefficients),space(u,2)).').'