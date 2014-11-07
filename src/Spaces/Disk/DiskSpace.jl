

include("JacobiSquareSpace.jl")

export Disk

##TODO: make argument 
immutable Disk <: BivariateDomain
    radius::Float64
    center::(Float64,Float64)
end

Disk(r)=Disk(r,(0.,0.))
Disk()=Disk(1.)

#canonical is rectangle [r,0]x[-π,π]
# we assume radius and centre are zero for now
fromcanonical(D::Disk,x,t)=x*cos(t),x*sin(t)
tocanonical(D::Disk,x,y)=sqrt(x^2+y^2),atan2(y,x)


# function points(d::Disk,n,m,k)
#     ptsx=0.5*(1-gaussjacobi(n,1.,0.)[1])
#     ptst=points(PeriodicInterval(),m)
#     
#     Float64[fromcanonical(d,x,t)[k] for x in ptsx, t in ptst]
# end


∂(d::Disk)=Circle(Complex(d.center...),d.radius)


# Kind== 0 => Legendre, K==1=>Chebyshev, K==2=>JacobiSquare
immutable DiskSpace{K,JS<:IntervalDomainSpace,S<:PeriodicDomainSpace} <: AbstractProductSpace{JS,S}
    domain::Disk
    spacet::S
end


DiskSpace{SS}(K::Integer,D::Disk,S::SS)=DiskSpace{K,K==2?JacobiSquareSpace:JacobiWeightSpace{JacobiSpace},SS}(D,S)
DiskSpace{SS}(D::Disk,S::SS)=DiskSpace{2,JacobiSquareSpace,SS}(D,S)
DiskSpace(K::Integer,D::Disk)=DiskSpace(K,D,LaurentSpace())
DiskSpace(K::Integer)=DiskSpace(K,Disk())

#TODO: Change to Fourier
DiskSpace(D::Disk)=DiskSpace(D,LaurentSpace())

domain(d::DiskSpace)=d.domain
function space(D::DiskSpace,k::Integer)
    @assert k==2
    D.spacet
end

Base.getindex(D::DiskSpace,k::Integer)=space(D,k)

Space(D::Disk)=DiskSpace(D)





columnspace{SS}(D::DiskSpace{0,SS},k)=(m=1.div(k,2);JacobiWeightSpace(0.,m,JacobiSpace(2m+1,0.,Interval(D.domain.radius,0.))))
columnspace{SS}(D::DiskSpace{1,SS},k)=(m=1.div(k,2);JacobiWeightSpace(0.,m,JacobiSpace(2m+0.5,-0.5,Interval(D.domain.radius,0.))))
columnspace{SS}(D::DiskSpace{2,SS},k)=(m=div(k,2);JacobiSquareSpace(m,Interval(D.domain.radius,0.)))

#transform(S::DiskSpace,V::Matrix)=transform([columnspace(S,k) for k=1:size(V,2)],S.spacet,V)


diskspacetype{K}(D::DiskSpace{K})=K


function Base.real{JS,D<:DiskSpace}(f::ProductFun{JS,LaurentSpace,D})
    cfs=f.coefficients
    n=length(cfs)

    ret=Array(Fun{JS,Float64},iseven(n)?n+1:n)
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

    ProductFun(ret,DiskSpace{diskspacetype(space(f)),JS,FourierSpace}(space(f).domain,FourierSpace()))
end
#Base.imag{S,T}(u::ProductFun{S,LarentSpace,T})=real(TensorFun(imag(u.coefficients),space(u,2)).').'+imag(TensorFun(real(u.coefficients),space(u,2)).').'




## Operators

function lap(S::Disk)
    D=Derivative()
    r=Fun(identity,[S.radius,0.])
    PDEOperator(((D^2+(1./r)*D)⊗I+(1./r).^2⊗D^2).ops,S)
end

neumann(S::Disk)=PDEOperator((lneumann()⊗I).ops,S)
dirichlet(S::Disk)=PDEOperator((ldirichlet()⊗I).ops,S)
diffbcs(S::Disk,k::Integer)=PDEOperator((ldiffbc(k)⊗I).ops,S)


