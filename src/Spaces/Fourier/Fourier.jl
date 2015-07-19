

export Fourier,Taylor,Hardy,CosSpace,SinSpace,Laurent

for T in (:CosSpace,:SinSpace)
    @eval begin
        immutable $T <: RealUnivariateSpace
            domain::PeriodicDomain
        end
        $T()=$T(PeriodicInterval())
        spacescompatible(a::$T,b::$T)=domainscompatible(a,b)
        canonicalspace(S::$T)=Fourier(domain(S))
    end
end

# s == true means analytic inside, taylor series
# s == false means anlytic outside and decaying at infinity
immutable Hardy{s} <: UnivariateSpace{ComplexBasis}
    domain::PeriodicDomain
    Hardy(d)=new(d)
    Hardy()=new(Circle())
end

canonicalspace(S::Hardy)=S

spacescompatible{s}(a::Hardy{s},b::Hardy{s})=domainscompatible(a,b)


typealias Taylor Hardy{true}

plan_transform(::Taylor,x::Vector)=wrap_fft_plan(plan_fft(x))
plan_itransform(::Taylor,x::Vector)=wrap_fft_plan(plan_ifft(x))
transform(::Taylor,vals::Vector,plan)=alternatesign!(plan(vals)/length(vals))
itransform(::Taylor,cfs::Vector,plan)=plan(alternatesign!(cfs))*length(cfs)

plan_transform(::Hardy{false},x::Vector)=wrap_fft_plan(plan_fft(x))
plan_itransform(::Hardy{false},x::Vector)=wrap_fft_plan(plan_ifft(x))
transform(::Hardy{false},vals::Vector,plan)=-alternatesign!(flipdim(plan(vals),1)/length(vals))
itransform(::Hardy{false},cfs::Vector,plan)=plan(flipdim(alternatesign!(-cfs),1))*length(cfs)

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
plan_transform(::CosSpace,x::Vector)=plan_chebyshevtransform(x;kind=2)
plan_itransform(::CosSpace,x::Vector)=plan_ichebyshevtransform(x;kind=2)
transform(::CosSpace,vals,plan)=chebyshevtransform(vals,plan;kind=2)
itransform(::CosSpace,cfs,plan)=ichebyshevtransform(cfs,plan;kind=2)
evaluate(f::Fun{CosSpace},t)=clenshaw(f.coefficients,cos(tocanonical(f,t)))


points(sp::SinSpace,n)=points(domain(sp),2n+2)[n+3:2n+2]
plan_transform{T<:FFTW.fftwNumber}(::SinSpace,x::Vector{T})=wrap_fft_plan(FFTW.plan_r2r(x,FFTW.RODFT00))
plan_itransform{T<:FFTW.fftwNumber}(::SinSpace,x::Vector{T})=wrap_fft_plan(FFTW.plan_r2r(x,FFTW.RODFT00))
transform(::SinSpace,vals,plan)=plan(vals)/(length(vals)+1)
itransform(::SinSpace,cfs,plan)=plan(cfs)/2
evaluate(f::Fun{SinSpace},t)=sineshaw(f.coefficients,tocanonical(f,t))



## Laurent space

typealias Laurent SumSpace{@compat(Tuple{Hardy{true},Hardy{false}}),ComplexBasis,1}

Laurent()=Laurent(PeriodicInterval())
Laurent{T<:Number}(d::Vector{T}) = Laurent(convert(PeriodicDomain,d))


points(sp::Laurent,n)=points(domain(sp),n)
plan_transform(::Laurent,x::Vector)=plan_svfft(x)
plan_itransform(::Laurent,x::Vector)=plan_isvfft(x)
transform(::Laurent,vals,plan)=svfft(vals,plan)
itransform(::Laurent,cfs,plan)=isvfft(cfs,plan)


## Fourier space

typealias Fourier SumSpace{@compat(Tuple{CosSpace,SinSpace}),RealBasis,1}
Fourier()=Fourier(PeriodicInterval())
Fourier{T<:Number}(d::Vector{T}) = Fourier(PeriodicInterval(d))


points(sp::Fourier,n)=points(domain(sp),n)
plan_transform{T<:FFTW.fftwNumber}(::Fourier,x::Vector{T}) = wrap_fft_plan(FFTW.plan_r2r(x, FFTW.R2HC))
plan_itransform{T<:FFTW.fftwNumber}(::Fourier,x::Vector{T}) = wrap_fft_plan(FFTW.plan_r2r(x, FFTW.HC2R))

function transform{T<:Number}(::Fourier,vals::Vector{T},plan)
    n=length(vals)
    cfs=2plan(vals)/n
    cfs[1]/=2
    if iseven(n)
        cfs[div(n,2)+1]/=2
    end

    fouriermodalt!(cfs)

    ret=Array(T,n)
    if iseven(n)
        ret[1:2:end]=cfs[1:div(n,2)]
        ret[2:2:end]=cfs[end:-1:div(n,2)+1]
    else
        ret[1:2:end]=cfs[1:div(n+1,2)]
        ret[2:2:end]=cfs[end:-1:div(n+3,2)]
    end
    ret
end

function itransform{T<:Number}(::Fourier,a::Vector{T},plan)
    n=length(a)
    cfs=[a[1:2:end];
            flipdim(a[2:2:end],1)]
    fouriermodalt!(cfs)
    if iseven(n)
        cfs[div(n,2)+1]*=2
    end
    cfs[1]*=2
    plan(cfs)/2
end

function fouriermodalt!(cfs)
    n=length(cfs)
    if iseven(n)
        for k=2:2:div(n,2)+1
            cfs[k]*=-1
        end
    else
        for k=2:2:div(n+1,2)
            cfs[k]*=-1
        end
    end

    if mod(n,4)==0
        for k=div(n,2)+3:2:n
            cfs[k]*=-1
        end
    elseif mod(n,4)==2
        for k=div(n,2)+2:2:n
            cfs[k]*=-1
        end
    elseif mod(n,4)==1
        for k=div(n+3,2):2:n
            cfs[k]*=-1
        end
    else #mod(n,4)==3
        for k=div(n+5,2):2:n
            cfs[k]*=-1
        end
    end
    cfs
end


Space(d::PeriodicInterval)=Fourier(d)
Space(d::Circle)=Laurent(d)
canonicalspace(S::Union(Laurent,Fourier))=isa(domain(S),Circle)?Laurent(domain(S)):Fourier(domain(S))

union_rule(A::ConstantSpace,B::Union(CosSpace,Taylor))=B

## Ones and zeros

for sp in (:Fourier,:Laurent,:Taylor,:CosSpace)
    @eval begin
        Base.ones{T<:Number}(::Type{T},S::$sp)=Fun(ones(T,1),S)
        Base.ones(S::$sp)=Fun(ones(1),S)
    end
end

reverseorientation(f::Fun{Fourier})=Fun(alternatesign!(copy(f.coefficients)),Fourier(reverse(domain(f))))

include("calculus.jl")
include("specialfunctions.jl")
include("FourierOperators.jl")
include("LaurentOperators.jl")
include("LaurentDirichlet.jl")




