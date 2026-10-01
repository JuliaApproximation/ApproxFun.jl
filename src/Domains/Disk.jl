export Disk,𝔻

# T is Real or Complex
# TT is (Real,Real) or Complex
immutable Disk{T,TT} <: BivariateDomain{T}
    center::TT
    radius::Float64
end

Disk(c::Complex,r)=Disk{typeof(c),typeof(c)}(c,r)
Disk{T<:Real}(c::Tuple{T,T},r)=Disk{T,typeof(c)}(c,r)

Disk()=Disk((0.,0.),1.)
Disk(::AnyDomain)=Disk(NaN,(NaN,NaN))

const 𝔻=Disk()

isambiguous(d::Disk)=isnan(d.radius) && all(isnan,d.center)


#canonical is rectangle [r,0]x[-π,π]
# we assume radius and centre are zero for now
fromcanonical{T<:Real}(D::Disk{T},x,t)=D.radius*x*cos(t)+D.center[1],D.radius*x*sin(t)+D.center[2]
tocanonical{T<:Real}(D::Disk{T},x,y)=sqrt((x-D.center[1])^2+(y-D.center[2])^2)/D.radius,atan2(y-D.center[2],x-D.center[1])
checkpoints(d::Disk)=[fromcanonical(d,(.1,.2243));fromcanonical(d,(-.212423,-.3))]

# function points(d::Disk,n,m,k)
#     ptsx=0.5*(1-gaussjacobi(n,1.,0.)[1])
#     ptst=points(PeriodicInterval(),m)
#
#     Float64[fromcanonical(d,x,t)[k] for x in ptsx, t in ptst]
# end


∂(d::Disk)=Circle(Complex(d.center...),d.radius)
