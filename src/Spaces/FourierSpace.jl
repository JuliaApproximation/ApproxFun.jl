

export FourierSpace,TaylorSpace,PoleSpace,CosSpace,SinSpace

for T in (:FourierSpace,:TaylorSpace,:PoleSpace,:CosSpace,:SinSpace)
    @eval begin
        immutable $T <: PeriodicDomainSpace
            domain::Union(PeriodicDomain,AnyDomain)
        end
        ==(a::($T),b::($T))= a.domain==b.domain
    end
end




transform(::TaylorSpace,vals::Vector)=alternatesign!(fft(vals)/length(vals))
itransform(::TaylorSpace,vals::Vector)=alternatesign!(ifft(vals)*length(vals))

function evaluate{T}(f::IFun{T,TaylorSpace},z)
    d=domain(f)
    if typeof(d) <: Circle
        horner(f.coefficients,(z-d.center)/d.radius)
    else
        horner(f.coefficients,fromcanonical(Circle(),tocanonical(f,z)))
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