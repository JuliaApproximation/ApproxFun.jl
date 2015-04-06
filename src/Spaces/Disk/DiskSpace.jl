

include("JacobiSquare.jl")

export Disk

##TODO: make argument
immutable Disk <: BivariateDomain{Float64}
    radius::Float64
    center::(Float64,Float64)
end

Disk(r)=Disk(r,(0.,0.))
Disk()=Disk(1.)

#canonical is rectangle [r,0]x[-π,π]
# we assume radius and centre are zero for now
fromcanonical(D::Disk,x,t)=x*cos(t),x*sin(t)
tocanonical(D::Disk,x,y)=sqrt(x^2+y^2),atan2(y,x)
checkpoints(d::Disk)=[fromcanonical(d,(.1,.2243));fromcanonical(d,(-.212423,-.3))]

# function points(d::Disk,n,m,k)
#     ptsx=0.5*(1-gaussjacobi(n,1.,0.)[1])
#     ptst=points(PeriodicInterval(),m)
#
#     Float64[fromcanonical(d,x,t)[k] for x in ptsx, t in ptst]
# end


∂(d::Disk)=Circle(Complex(d.center...),d.radius)


# Kind== 0 => Legendre, K==1=>Chebyshev, K==2=>JacobiSquare
immutable DiskSpace{m,a,b,JS<:IntervalSpace,S<:PeriodicSpace} <: AbstractProductSpace{JS,S,Complex128,Disk}
    domain::Disk
    spacet::S
end


DiskSpace(m,a,b,D::Disk,S::PeriodicSpace)=DiskSpace{m,a,b,JacobiSquare,typeof(S)}(D,S)
DiskSpace(D::Disk,S::PeriodicSpace)=DiskSpace(0,0,0,D,S)
DiskSpace(D::Disk)=DiskSpace(D,Laurent())

spacescompatible{m,a,b,JS,S}(A::DiskSpace{m,a,b,JS,S},B::DiskSpace{m,a,b,JS,S})=true

coefficient_type{T<:Complex}(::DiskSpace,::Type{T})=T
coefficient_type{T<:Real}(::DiskSpace,::Type{T})=Complex{T}

domain(d::DiskSpace)=d.domain
function space(D::DiskSpace,k::Integer)
    @assert k==2
    D.spacet
end

Base.getindex(D::DiskSpace,k::Integer)=space(D,k)

Space(D::Disk)=DiskSpace(D)


columnspace{M,a,b,SS}(D::DiskSpace{M,a,b,SS},k)=(m=div(k,2);JacobiSquare(M+m,a+m,b,Interval(D.domain.radius,0.)))

#transform(S::DiskSpace,V::Matrix)=transform([columnspace(S,k) for k=1:size(V,2)],S.spacet,V)


function Base.real{JS}(f::ProductFun{JS,Laurent,DiskSpace{0,0,0,JS,Laurent}})
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

    ProductFun(ret,DiskSpace{0,0,0,JS,Fourier}(space(f).domain,Fourier()))
end
#Base.imag{S,T}(u::ProductFun{S,Larent,T})=real(TensorFun(imag(u.coefficients),space(u,2)).').'+imag(TensorFun(real(u.coefficients),space(u,2)).').'



## Conversion
# These are placeholders for future

conversion_rule{m,a,b,m2,a2,b2,JS,FS}(A::DiskSpace{m,a,b,JS,FS},B::DiskSpace{m2,a2,b2,JS,FS})=DiskSpace(max(m,m2),min(a,a2),min(b,b2),A.domain,B.spacet)


## Operators

isfunctional{DS<:DiskSpace}(::Dirichlet{DS},k)=k==1
dekron{DS<:DiskSpace}(::Dirichlet{DS},k)=k==1?ldirichlet():ConstantOperator(1.0)


Base.length{DS<:DiskSpace}(::Laplacian{DS})=2
function dekron{DS<:DiskSpace}(L::Laplacian{DS},k,::Colon)
    if k==1
        r=Fun(identity,[domain(L).radius,0.])
        D=Derivative()
        [(D^2+(1./r)*D),Multiplication((1./r).^2)]
    elseif k==2
        D=Derivative(domainspace(L)[2])
        [ConstantOperator(1.0),D^2]
    end
end

dekron(L,k::Integer,j::Integer)=dekron(L,k,:)[j]

lap(d::Disk)=Laplacian(Space(d))
dirichlet(d::Disk)=Dirichlet(Space(d))



