

export Fourier,Taylor,Hardy,CosSpace,SinSpace,Laurent

for T in (:CosSpace,:SinSpace)
    @eval begin
        immutable $T <: PeriodicSpace{Float64}
            domain::Union(PeriodicDomain,AnyDomain)
        end
        $T()=$T(PeriodicInterval())
    end
end

# s == true means analytic inside, taylor series
# s == false means anlytic outside and decaying at infinity
immutable Hardy{s} <: PeriodicSpace{Complex{Float64}}
    domain::Union(PeriodicDomain,AnyDomain)
    Hardy(d)=new(d)
    Hardy()=new(Circle())
end



typealias Taylor Hardy{true}

transform(::Taylor,vals::Vector)=alternatesign!(fft(vals)/length(vals))
itransform(::Taylor,cfs::Vector)=ifft(alternatesign!(cfs))*length(cfs)

transform(::Hardy{false},vals::Vector)=-alternatesign!(flipud(fft(vals))/length(vals))
itransform(::Hardy{false},cfs::Vector)=ifft(flipud(alternatesign!(-cfs)))*length(cfs)

function evaluate(f::Fun{Taylor},z)
    d=domain(f)
    if isa(d,Circle)
        horner(f.coefficients,(z-d.center)/d.radius)
    else
        horner(f.coefficients,fromcanonical(Circle(),tocanonical(f,z)))
    end
end

function evaluate(f::Fun{Hardy{false}},z)
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
function horner{T}(v::Vector{T},z::Number)
    ret = zero(T)
    ei = z*one(T)
    
    p = one(T)
    for vk in v
        ret += vk*p
        p .*= ei
    end
    
    ret
end

function horner{T}(v::Vector{T},z::Vector)
    ret = zeros(T,length(z))
    ei = z*one(T)
    
    p = ones(T,length(z))
    for vk in v
        ret += vk*p
        p .*= ei
    end
    
    ret
end

## Cos and Sin space

points(sp::CosSpace,n)=points(domain(sp),2n-2)[1:n]
transform(::CosSpace,vals)=chebyshevtransform(vals)
itransform(::CosSpace,cfs)=ichebyshevtransform(cfs)
evaluate(f::Fun{CosSpace},t)=clenshaw(f.coefficients,cos(tocanonical(f,t)))


points(sp::SinSpace,n)=fromcanonical(domain(sp),(π*[1:n])/(n+1))
transform(::SinSpace,vals)=FFTW.r2r(vals,FFTW.RODFT00)/(length(vals)+1)
itransform(::SinSpace,cfs)=FFTW.r2r(cfs,FFTW.RODFT00)/2
evaluate{T,N<:Number}(f::Fun{SinSpace,T},t::N)=sum(promote_type(T,N)[f.coefficients[k]*sin(k*tocanonical(f,t)) for k=1:length(f)])



## Laurent space

typealias Laurent PeriodicSumSpace{Hardy{true},Hardy{false},Complex{Float64}}
Laurent()=Laurent(PeriodicInterval())
Laurent{T<:Number}(d::Vector{T}) = Laurent(PeriodicInterval(d))

Space(d::PeriodicDomain)=Fourier(d)
Space(d::Circle)=Laurent(d)
canonicalspace(S::PeriodicSpace)=isa(domain(S),Circle)?Laurent(domain(S)):Fourier(domain(S))


points(sp::Laurent,n)=points(domain(sp),n)
transform(::Laurent,vals)=svfft(vals)|>interlace
itransform(::Laurent,cfs)=isvfft(deinterlace(cfs))

## Ones and zeros


Base.ones{T<:Number}(::Type{T},S::Laurent)=Fun(ones(T,1),S)


## Fourier space

typealias Fourier PeriodicSumSpace{CosSpace,SinSpace,Float64}
Fourier()=Fourier(PeriodicInterval())
Fourier{T<:Number}(d::Vector{T}) = Fourier(PeriodicInterval(d))


#domain(S) may be any domain
for sp in (:Fourier,:Laurent,:(Hardy{true}),:CosSpace)
    @eval begin
        Base.ones{T<:Number}(::Type{T},S::$sp)=Fun(ones(T,1),S)
        Base.ones(S::$sp)=Fun(ones(1),S)        
    end
end

points(sp::Fourier,n)=points(domain(sp),n)


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

function transform{T<:Number}(::Fourier,vals::Vector{T})
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


function itransform{T<:Number}(::Fourier,a::Vector{T})
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
include("LaurentOperators.jl")