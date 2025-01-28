

export Fourier,Taylor,Hardy,CosSpace,SinSpace,Laurent

for T in (:CosSpace,:SinSpace)
    @eval begin
        doc"""
            `CosSpace` is the basis `[1,cos θ,cos 2θ,...]`
            `SinSpace` is the basis `[sin θ,sin 2θ,...]`
        """
        immutable $T{D<:Domain} <: RealUnivariateSpace{D}
            domain::D
            $T(d::Domain) = new(D(d))
            $T(d::D) = new(d)
        end
        $T(d::Domain) = $T{typeof(d)}(d)
        $T() = $T(PeriodicInterval())
        spacescompatible(a::$T,b::$T) = domainscompatible(a,b)
        hasfasttransform(::$T) = true
        canonicalspace(S::$T) = Fourier(domain(S))
        setdomain(S::$T,d::Domain) = $T(d)
    end
end
# s == true means analytic inside, taylor series
# s == false means anlytic outside and decaying at infinity

doc"""
    `Hardy{true}` is the basis `[1,z,z^2,...]`
    `Hardy{false}` is the basis `[1/z,1/z^2,...]`
"""
immutable Hardy{s,D<:Domain} <: UnivariateSpace{ComplexBasis,D}
    domain::D
    Hardy(d)=new(d)
    Hardy()=new(D())
end


Base.promote_rule{T<:Number,S<:Union{Hardy{true},CosSpace},V}(::Type{Fun{S,V}},::Type{T}) =
    Fun{S,promote_type(V,T)}
Base.promote_rule{T<:Number,S<:Union{Hardy{true},CosSpace}}(::Type{Fun{S}},::Type{T}) =
    Fun{S,T}

@compat (H::Type{Hardy{s}}){s}(d::Domain) = Hardy{s,typeof(d)}(d)
@compat (H::Type{Hardy{s}}){s}() = Hardy{s}(Circle())

canonicalspace(S::Hardy) = S
setdomain{s}(S::Hardy{s},d::Domain) = Hardy{s}(d)


spacescompatible{s}(a::Hardy{s},b::Hardy{s}) = domainscompatible(a,b)
hasfasttransform(::Hardy) = true

# The <: Domain is crucial for matching Basecall overrides
typealias Taylor{D<:Domain} Hardy{true,D}

plan_transform(::Taylor,x::Vector) = plan_fft(x)
plan_itransform(::Taylor,x::Vector) = plan_ifft(x)
transform(::Taylor,vals::Vector,plan) = alternatesign!(plan*vals/length(vals))
itransform(::Taylor,cfs::Vector,plan) = plan*alternatesign!(cfs)*length(cfs)

plan_transform(::Hardy{false},x::Vector) = plan_fft(x)
plan_itransform(::Hardy{false},x::Vector) = plan_ifft(x)
transform(::Hardy{false},vals::Vector,plan) = -alternatesign!(flipdim(plan*vals,1)/length(vals))
itransform(::Hardy{false},cfs::Vector,plan) = plan*flipdim(alternatesign!(-cfs),1)*length(cfs)

evaluate{D<:Domain}(f::AbstractVector,S::Taylor{D},z) = horner(f,fromcanonical(Circle(),tocanonical(S,z)))
function evaluate{D<:Circle}(f::AbstractVector,S::Taylor{D},z)
    z=mappoint(S,𝕌,z)
    d=domain(S)
    horner(f,z)
end

function evaluate{D<:Domain}(f::AbstractVector,S::Hardy{false,D},z)
    z=mappoint(S,𝕌,z)
    z=1./z
    z.*horner(f,z)
end
function evaluate{D<:Circle}(f::AbstractVector,S::Hardy{false,D},z)
    z=mappoint(S,𝕌,z)
    z=1./z
    z.*horner(f,z)
end


##TODO: fast routine

function horner(c::AbstractVector,kr::Range{Int64},x)
    T = promote_type(eltype(c),eltype(x))
    if isempty(c)
        return zero(x)
    end

    ret = zero(T)
    @inbounds for k in reverse(kr)
        ret = muladd(x,ret,c[k])
    end

    ret
end

function horner(c::AbstractVector,kr::Range{Int64},x::AbstractVector)
    n,T = length(x),promote_type(eltype(c),eltype(x))
    if isempty(c)
        return zero(x)
    end

    ret = zeros(T,n)
    @inbounds for k in reverse(kr)
        ck = c[k]
        @simd for i = 1:n
            ret[i] = muladd(x[i],ret[i],ck)
        end
    end

    ret
end

horner(c::AbstractVector,x) = horner(c,1:length(c),x)
horner(c::AbstractVector,x::AbstractArray) = horner(c,1:length(c),x)
horner(c::AbstractVector,kr::Range{Int64},x::AbstractArray) = reshape(horner(c,kr,vec(x)),size(x))

## Cos and Sin space

points(sp::CosSpace,n)=points(domain(sp),2n-2)[1:n]
plan_transform(::CosSpace,x::Vector)=plan_chebyshevtransform(x;kind=2)
plan_itransform(::CosSpace,x::Vector)=plan_ichebyshevtransform(x;kind=2)
transform(::CosSpace,vals,plan)=chebyshevtransform(vals,plan;kind=2)
itransform(::CosSpace,cfs,plan)=ichebyshevtransform(cfs,plan;kind=2)
evaluate(f::Vector,S::CosSpace,t)=clenshaw(Chebyshev(),f,cos(tocanonical(S,t)))


points(sp::SinSpace,n)=points(domain(sp),2n+2)[n+3:2n+2]
plan_transform{T<:FFTW.fftwNumber}(::SinSpace,x::Vector{T})=FFTW.plan_r2r(x,FFTW.RODFT00)
plan_itransform{T<:FFTW.fftwNumber}(::SinSpace,x::Vector{T})=FFTW.plan_r2r(x,FFTW.RODFT00)

plan_transform{D}(::SinSpace{D},x::Vector)=error("transform for Fourier only implemented for fftwNumbers")
plan_itransform{D}(::SinSpace{D},x::Vector)=error("transform for Fourier only implemented for fftwNumbers")



transform(::SinSpace,vals,plan)=plan*vals/(length(vals)+1)
itransform(::SinSpace,cfs,plan)=plan*cfs/2
evaluate(f::AbstractVector,S::SinSpace,t)=sineshaw(f,tocanonical(S,t))



## Laurent space

typealias Laurent{DD} SumSpace{Tuple{Hardy{true,DD},Hardy{false,DD}},ComplexBasis,DD,1}


plan_transform{DD}(::Laurent{DD},x::Vector)=plan_svfft(x)
plan_itransform{DD}(::Laurent{DD},x::Vector)=plan_isvfft(x)
transform{DD}(::Laurent{DD},vals,plan)=svfft(vals,plan)
itransform{DD}(::Laurent{DD},cfs,plan)=isvfft(cfs,plan)


function evaluate{DD}(f::AbstractVector,S::Laurent{DD},z)
    z = mappoint(domain(S),Circle(),z)
    invz = 1./z
    horner(f,1:2:length(f),z) + horner(f,2:2:length(f),invz).*invz
end


function Base.conj{DD}(f::Fun{Laurent{DD}})
    cfs=Array(eltype(f),iseven(ncoefficients(f))?ncoefficients(f)+1:ncoefficients(f))
    cfs[1]=conj(f.coefficients[1])
    for k=2:2:ncoefficients(f)-1
        cfs[k]=conj(f.coefficients[k+1])
    end
    for k=3:2:ncoefficients(f)
        cfs[k]=conj(f.coefficients[k-1])
    end
    Fun(cfs,space(f))
end

## Fourier space

typealias Fourier{DD} SumSpace{Tuple{CosSpace{DD},SinSpace{DD}},RealBasis,DD,1}

for TYP in (:Laurent,:Fourier)
    @eval begin
        @compat (::Type{$TYP})(d::Domain) = $TYP{typeof(d)}(d)
        @compat (::Type{$TYP})() = $TYP(PeriodicInterval())
        @compat (::Type{$TYP}){T<:Number}(d::Vector{T}) = $TYP(PeriodicDomain(d))

        hasfasttransform{D}(::$TYP{D}) = true
    end
end

for T in (:CosSpace,:SinSpace)
    @eval begin
        # override default as canonicalspace must be implemented
        maxspace{D}(::$T,::Fourier{D}) = NoSpace()
        maxspace{D}(::Fourier{D},::$T) = NoSpace()
    end
end

points{D}(sp::Fourier{D},n)=points(domain(sp),n)
plan_transform{T<:FFTW.fftwNumber,D}(::Fourier{D},x::Vector{T}) = FFTW.plan_r2r(x, FFTW.R2HC)
plan_itransform{T<:FFTW.fftwNumber,D}(::Fourier{D},x::Vector{T}) = FFTW.plan_r2r(x, FFTW.HC2R)

plan_transform{D}(::Fourier{D},x::Vector)=error("transform for Fourier only implemented for fftwNumbers")
plan_itransform{D}(::Fourier{D},x::Vector)=error("transform for Fourier only implemented for fftwNumbers")

function transform{T<:Number,D}(::Fourier{D},vals::Vector{T},plan)
    n=length(vals)
    cfs=2plan*vals/n
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

function itransform{T<:Number,D}(::Fourier{D},a::Vector{T},plan)
    n=length(a)
    cfs=[a[1:2:end];
            flipdim(a[2:2:end],1)]
    fouriermodalt!(cfs)
    if iseven(n)
        cfs[div(n,2)+1]*=2
    end
    cfs[1]*=2
    plan*cfs/2
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



canonicalspace{DD<:PeriodicInterval}(S::Laurent{DD})=Fourier(domain(S))
canonicalspace{DD<:Circle}(S::Fourier{DD})=Laurent(domain(S))
canonicalspace{DD<:PeriodicLine}(S::Laurent{DD})=S

for TYP in (:CosSpace,:Taylor)
    @eval union_rule(A::ConstantSpace,B::$TYP)=B
end

## Ones and zeros

for sp in (:Fourier,:CosSpace,:Laurent,:Taylor)
    @eval begin
        Base.ones{T<:Number,D}(::Type{T},S::$sp{D})=Fun(ones(T,1),S)
        Base.ones{D}(S::$sp{D})=Fun(ones(1),S)
    end
end


function identity_fun{DD<:Circle}(S::Taylor{DD})
    d=domain(S)
    if d.orientation
        Fun([d.center,d.radius],S)
    else
        error("Cannot create identity on $S")
    end
end


identity_fun{DD<:Circle}(S::Fourier{DD}) = Fun(identity_fun(Laurent(domain(S))),S)


reverseorientation{D}(f::Fun{Fourier{D}})=Fun(alternatesign!(copy(f.coefficients)),Fourier(reverse(domain(f))))
function reverseorientation{D}(f::Fun{Laurent{D}})
    # exp(im*k*x) -> exp(-im*k*x), or equivalentaly z -> 1/z
    n=ncoefficients(f)
    ret=Array(eltype(f),iseven(n)?n+1:n)  # since z -> 1/z we get one more coefficient
    ret[1]=f.coefficients[1]
    for k=2:2:length(ret)-1
        ret[k+1]=f.coefficients[k]
    end
    for k=2:2:n-1
        ret[k]=f.coefficients[k+1]
    end
    iseven(n) && (ret[n] = 0)

    Fun(ret,Laurent(reverse(domain(f))))
end

include("calculus.jl")
include("specialfunctions.jl")
include("FourierOperators.jl")
include("LaurentOperators.jl")
include("LaurentDirichlet.jl")
