export Disk,𝔻

# T is Real or Complex
# TT is (Real,Real) or Complex
"""
    Disk(c,r)

represents the disk centred at `c` with radius `r`.
"""
struct Disk{T} <: Domain{T}
    center::T
    radius::Float64
end

Disk(c,r) = Disk{typeof(c)}(c,r)
Disk(c::Real,r) = Disk(complex(c),r)
Disk(c::Integer,r) = Disk(float(c),r)
Disk(c::Complex{<:Integer},r) = Disk(float(c),r)

Disk() = Disk(Vec(0.,0.),1.)
Disk(::AnyDomain) = Disk(NaN,(NaN,NaN))

const 𝔻 = Disk()

isambiguous(d::Disk) = isnan(d.radius) && all(isnan,d.center)


#canonical is rectangle [r,0]x[-π,π]
# we assume radius and centre are zero for now
fromcanonical(D::Disk{T},x,t) where {T<:Real} =
    Vec(D.radius*x*cos(t)+D.center[1],D.radius*x*sin(t)+D.center[2])
tocanonical(D::Disk{T},x,y) where {T<:Real} =
    Vec(sqrt((x-D.center[1])^2+(y-D.center[2])^2)/D.radius,
        atan(y-D.center[2],x-D.center[1]))
checkpoints(d::Disk) = [fromcanonical(d,(.1,.2243));fromcanonical(d,(-.212423,-.3))]

# function points(d::Disk,n,m,k)
#     ptsx=0.5*(1-gaussjacobi(n,1.,0.)[1])
#     ptst=points(PeriodicSegment(),m)
#
#     Float64[fromcanonical(d,x,t)[k] for x in ptsx, t in ptst]
# end


boundary(d::Disk) = Circle(d.center,d.radius)
