

export FourierSpace,TaylorSpace,HardySpace,CosSpace,SinSpace,LaurentSpace

for T in (:FourierSpace,:CosSpace,:SinSpace)
    @eval begin
        immutable $T <: PeriodicDomainSpace
            domain::Union(PeriodicDomain,AnyDomain)
        end
        ==(a::($T),b::($T))= a.domain==b.domain
    end
end

# s == true means analytic inside, taylor series
# s == false means anlytic outside and decaying at infinity
immutable HardySpace{s} <: PeriodicDomainSpace
    domain::Union(PeriodicDomain,AnyDomain)
end

=={s}(a::HardySpace{s},b::HardySpace{s})= a.domain==b.domain

typealias TaylorSpace HardySpace{true}
typealias PoleSpace HardySpace{false}

transform(::TaylorSpace,vals::Vector)=alternatesign!(fft(vals)/length(vals))
itransform(::TaylorSpace,cfs::Vector)=ifft(alternatesign!(cfs))*length(cfs)

transform(::PoleSpace,vals::Vector)=-alternatesign!(flipud(fft(vals))/length(vals))
itransform(::PoleSpace,cfs::Vector)=ifft(flipud(alternatesign!(-cfs)))*length(cfs)

function evaluate{T}(f::IFun{T,TaylorSpace},z)
    d=domain(f)
    if typeof(d) <: Circle
        horner(f.coefficients,(z-d.center)/d.radius)
    else
        horner(f.coefficients,fromcanonical(Circle(),tocanonical(f,z)))
    end
end

function evaluate{T}(f::IFun{T,PoleSpace},z)
    d=domain(f)
    if typeof(d) <: Circle
        z=(z-d.center)/d.radius
        z=1./z
        z.*horner(f.coefficients,z)
    else
        z=fromcanonical(Circle(),tocanonical(f,z))
        z=1./z
        z.*horner(f.coefficients,z)
    end
end


##TODO: fast routine
function horner{T}(v::Vector{T},z)
    ret = zero(T)
    ei = z*one(T)
    
    p = one(T)
    for vk in v
        ret += vk*p
        p *= ei
    end
    
    ret
end

points(sp::CosSpace,n)=points(domain(sp),2n-2)[1:n]
transform(::CosSpace,vals)=chebyshevtransform(vals)
itransform(::CosSpace,vals)=ichebyshevtransform(vals)

evaluate{T}(f::IFun{T,CosSpace},t)=clenshaw(f.coefficients,cos(tocanonical(f,t)))



## Laurent space

typealias LaurentSpace SumSpace{HardySpace{true},HardySpace{false}}
LaurentSpace(d::PeriodicDomain)=SumSpace((HardySpace{true}(d),HardySpace{false}(d)))

points(sp::LaurentSpace,n)=points(domain(sp),n)
transform(::LaurentSpace,vals)=svfft(vals)|>interlace
itransform(::LaurentSpace,cfs)=isvfft(deinterlace(cfs))

