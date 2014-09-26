export Disk

immutable Disk <: BivariateDomain
    radius::Float64
    center::(Float64,Float64)
end
Disk()=Disk(0.,(0.,0.))

#canonical is rectangle [1,0]x[-π,π]
# we assume radius and centre are zero for now
fromcanonical(D::Disk,x,t)=x*cos(t),x*sin(t)
tocanonical(D::Disk,x,y)=sqrt(x^2+y^2),atan2(y,x)


function points(d::Disk,n,m,k)
    ptsx=0.5*(1-gaussjacobi(n,1.,0.)[1])
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

Base.getindex(D::DiskSpace,k::Integer)=space(D,k)

Space(D::Disk)=DiskSpace(D)

points(d::DiskSpace,n...)=points(domain(d),n...)


columnspace(D::DiskSpace,k)=(m=1.div(k,2);JacobiWeightSpace(0.,m,JacobiSpace(2m+1,0.,Interval(1.,0.))))

transform(S::DiskSpace,V::Matrix)=transform([columnspace(S,k) for k=1:size(V,2)],S.spacet,V)


function Base.real(f::ProductFun{JacobiWeightSpace{JacobiSpace},LaurentSpace,DiskSpace{LaurentSpace}})
    cfs=f.coefficients
    n=length(cfs)

    ret=Array(Fun{JacobiWeightSpace{JacobiSpace},Float64},iseven(n)?n+1:n)
    ret[1]=real(cfs[1])
    
    for k=2:2:n
        # exp(1im(k-1)/2*x)=cos((k-1)/2 x) +i sin((k-1)/2 x)
        ret[k]=imag(cfs[k])
        ret[k+1]=real(cfs[k])
    end        
    for k=3:2:n
        # exp(1im(k-1)/2*x)=cos((k-1)/2 x) +i sin((k-1)/2 x)
        ret[k]+=real(cfs[k])
        ret[k-1]-=imag(cfs[k])
    end

    ProductFun(ret,DiskSpace(space(f).domain,FourierSpace()))
end
#Base.imag{S,T}(u::ProductFun{S,LarentSpace,T})=real(TensorFun(imag(u.coefficients),space(u,2)).').'+imag(TensorFun(real(u.coefficients),space(u,2)).').'