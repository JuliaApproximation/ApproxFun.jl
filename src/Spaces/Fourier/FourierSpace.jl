

export FourierSpace,TaylorSpace,HardySpace,CosSpace,SinSpace,LaurentSpace

for T in (:CosSpace,:SinSpace)
    @eval begin
        immutable $T <: PeriodicDomainSpace{Float64}
            domain::Union(PeriodicDomain,AnyDomain)
        end
        ==(a::($T),b::($T))= a.domain==b.domain
    end
end

# s == true means analytic inside, taylor series
# s == false means anlytic outside and decaying at infinity
immutable HardySpace{s} <: PeriodicDomainSpace{Complex{Float64}}
    domain::Union(PeriodicDomain,AnyDomain)
end

=={s}(a::HardySpace{s},b::HardySpace{s})= a.domain==b.domain

typealias TaylorSpace HardySpace{true}
typealias PoleSpace HardySpace{false}

transform(::TaylorSpace,vals::Vector)=alternatesign!(fft(vals)/length(vals))
itransform(::TaylorSpace,cfs::Vector)=ifft(alternatesign!(cfs))*length(cfs)

transform(::PoleSpace,vals::Vector)=-alternatesign!(flipud(fft(vals))/length(vals))
itransform(::PoleSpace,cfs::Vector)=ifft(flipud(alternatesign!(-cfs)))*length(cfs)

function evaluate{T}(f::Fun{T,TaylorSpace},z)
    d=domain(f)
    if isa(d,Circle)
        horner(f.coefficients,(z-d.center)/d.radius)
    else
        horner(f.coefficients,fromcanonical(Circle(),tocanonical(f,z)))
    end
end

function evaluate{T}(f::Fun{T,PoleSpace},z)
    d=domain(f)
    if isa(d,Circle)
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

## Cos and Sin space

points(sp::CosSpace,n)=points(domain(sp),2n-2)[1:n]
transform(::CosSpace,vals)=chebyshevtransform(vals)
itransform(::CosSpace,cfs)=ichebyshevtransform(cfs)
evaluate{T}(f::Fun{T,CosSpace},t)=clenshaw(f.coefficients,cos(tocanonical(f,t)))


points(sp::SinSpace,n)=fromcanonical(domain(sp),(π*[1:n])/(n+1))
transform(::SinSpace,vals)=FFTW.r2r(vals,FFTW.RODFT00)/(length(vals)+1)
itransform(::SinSpace,cfs)=FFTW.r2r(cfs,FFTW.RODFT00)/2
evaluate{T}(f::Fun{T,SinSpace},t)=sum(T[f.coefficients[k]*sin(k*tocanonical(f,t)) for k=1:length(f)])



## Laurent space

typealias LaurentSpace PeriodicSumSpace{Complex{Float64},HardySpace{true},HardySpace{false}}
LaurentSpace(d::Union(PeriodicDomain,AnyDomain))=PeriodicSumSpace(HardySpace{true}(d),HardySpace{false}(d))

Space(d::PeriodicDomain)=LaurentSpace(d)
canonicalspace(S::PeriodicDomainSpace)=LaurentSpace(domain(S))


points(sp::LaurentSpace,n)=points(domain(sp),n)
transform(::LaurentSpace,vals)=svfft(vals)|>interlace
itransform(::LaurentSpace,cfs)=isvfft(deinterlace(cfs))

## Ones and zeros


for op in (:(Base.ones),:(Base.zeros))
    @eval ($op){T<:Number}(::Type{T},S::LaurentSpace)=Fun(($op)(T,1),S)
end


## Fourier space

typealias FourierSpace PeriodicSumSpace{Float64,CosSpace,SinSpace}
FourierSpace(d::Union(PeriodicDomain,AnyDomain))=PeriodicSumSpace((CosSpace(d),SinSpace(d)))

points(sp::FourierSpace,n)=points(domain(sp),n)


function fouriermodalt!(cfs)
    n=length(cfs)
    if iseven(n)
        for k=2:2:n/2+1
            cfs[k]*=-1
        end  
    else
        for k=2:2:(n+1)/2
            cfs[k]*=-1
        end     
    end    
    
    if mod(n,4)==0
        for k=n/2+3:2:n
            cfs[k]*=-1
        end
    elseif mod(n,4)==2
        for k=n/2+2:2:n
            cfs[k]*=-1
        end    
    elseif mod(n,4)==1
        for k=(n+3)/2:2:n
            cfs[k]*=-1
        end   
    else #mod(n,4)==3
        for k=(n+5)/2:2:n
            cfs[k]*=-1
        end      
    end
    cfs
end

function transform{T<:Number}(::FourierSpace,vals::Vector{T})
    n=length(vals)
    cfs=2FFTW.r2r(vals, FFTW.R2HC )/n
    cfs[1]/=2
    if iseven(n)
        cfs[n/2+1]/=2
    end

    fouriermodalt!(cfs)
        
    ret=Array(T,n)
    if iseven(n)
        ret[1:2:end]=cfs[1:n/2]
        ret[2:2:end]=cfs[end:-1:n/2+1]
    else
        ret[1:2:end]=cfs[1:(n+1)/2]
        ret[2:2:end]=cfs[end:-1:(n+3)/2]
    end
    ret    
end


function itransform{T<:Number}(::FourierSpace,a::Vector{T})
    n=length(a)
    cfs=[a[1:2:end],flipud(a[2:2:end])]
    fouriermodalt!(cfs)
    if iseven(n)
        cfs[n/2+1]*=2
    end        
    cfs[1]*=2
    FFTW.r2r(cfs, FFTW.HC2R )/2  
end









include("calculus.jl")
include("specialfunctions.jl")
include("FourierOperators.jl")