

export Fourier,Taylor,Hardy,CosSpace,SinSpace,Laurent

for T in (:CosSpace,:SinSpace)
    @eval begin
        struct $T{D<:Domain,R} <: Space{D,R}
            domain::D
            $T{D,R}(d::Domain) where {D,R} = new(D(d))
            $T{D,R}(d::D) where{D,R} = new(d)
        end
        $T(d::Domain) = $T{typeof(d),real(prectype(d))}(d)
        $T() = $T(PeriodicInterval())
        spacescompatible(a::$T,b::$T) = domainscompatible(a,b)
        hasfasttransform(::$T) = true
        canonicalspace(S::$T) = Fourier(domain(S))
        setdomain(S::$T,d::Domain) = $T(d)
    end
end

doc"""
`CosSpace()` is the space spanned by `[1,cos θ,cos 2θ,...]`
"""
CosSpace()

doc"""
`SinSpace()` is the space spanned by `[sin θ,sin 2θ,...]`
"""
SinSpace()

# s == true means analytic inside, taylor series
# s == false means anlytic outside and decaying at infinity
doc"""
`Hardy{false}()` is the space spanned by `[1/z,1/z^2,...]`.
`Hardy{true}()` is the space spanned by `[1,z,z^2,...]`.
"""
struct Hardy{s,D<:Domain,R} <: Space{D,R}
    domain::D
    Hardy{s,D,R}(d) where {s,D,R} = new{s,D,R}(d)
    Hardy{s,D,R}() where {s,D,R}  = new{s,D,R}(D())
end

# The <: Domain is crucial for matching Basecall overrides
doc"""
`Taylor()` is the space spanned by `[1,z,z^2,...]`.
This is a type alias for `Hardy{true}`.
"""
const Taylor{D<:Domain,R} = Hardy{true,D,R}


@containsconstants CosSpace
@containsconstants Taylor



Base.promote_rule{T<:Number,S<:Union{Hardy{true},CosSpace},V}(::Type{Fun{S,V}},::Type{T}) =
    Fun{S,promote_type(V,T)}
Base.promote_rule{T<:Number,S<:Union{Hardy{true},CosSpace}}(::Type{Fun{S}},::Type{T}) =
    Fun{S,T}

(H::Type{Hardy{s}}){s}(d::Domain) = Hardy{s,typeof(d),complex(prectype(d))}(d)
(H::Type{Hardy{s}}){s}() = Hardy{s}(Circle())

canonicalspace(S::Hardy) = S
setdomain{s}(S::Hardy{s},d::Domain) = Hardy{s}(d)


spacescompatible{s}(a::Hardy{s},b::Hardy{s}) = domainscompatible(a,b)
hasfasttransform(::Hardy) = true


for (Typ,Plfft!,Plfft,Pltr!,Pltr) in ((:TransformPlan,:plan_fft!,:plan_fft,:plan_transform!,:plan_transform),
                           (:ITransformPlan,:plan_ifft!,:plan_ifft,:plan_itransform!,:plan_itransform))
    @eval begin
        $Pltr!{T<:Complex}(sp::Hardy,x::AbstractVector{T}) = $Typ(sp,$Plfft!(x),Val{true})
        $Pltr!{T<:Real}(::Hardy,x::AbstractVector{T}) =
            error("In place variants not possible with real data.")

        $Pltr{T<:Complex}(sp::Hardy,x::AbstractVector{T}) = $Typ(sp,$Pltr!(sp,x),Val{false})
        function $Pltr{T}(sp::Hardy,x::AbstractVector{T})
            plan = $Pltr(sp,Array{Complex{T}}(length(x))) # we can reuse vector in itransform
            $Typ{T,typeof(sp),false,typeof(plan)}(sp,plan)
        end

        *{T<:Complex,HS<:Hardy}(P::$Typ{T,HS,false},vals::AbstractVector{T}) = P.plan*copy(vals)
        *{T,HS<:Hardy}(P::$Typ{T,HS,false},vals::AbstractVector{T}) = P.plan*Vector{Complex{T}}(vals)
    end
end


*{T,DD,RR}(P::TransformPlan{T,Hardy{true,DD,RR},true},vals::AbstractVector{T}) =
    scale!(one(T)/length(vals),P.plan*vals)
*{T,DD,RR}(P::ITransformPlan{T,Hardy{true,DD,RR},true},cfs::AbstractVector{T}) =
    scale!(length(cfs),P.plan*cfs)
*{T,DD,RR}(P::TransformPlan{T,Hardy{false,DD,RR},true},vals::AbstractVector{T}) =
    scale!(one(T)/length(vals),reverse!(P.plan*vals))
*{T,DD,RR}(P::ITransformPlan{T,Hardy{false,DD,RR},true},cfs::AbstractVector{T}) =
    scale!(length(cfs),P.plan*reverse!(cfs))


transform(sp::Hardy,vals::AbstractVector,plan) = plan*vals
itransform(sp::Hardy,vals::AbstractVector,plan) = plan*vals

evaluate{D<:Domain,R}(f::AbstractVector,S::Taylor{D,R},z) = horner(f,fromcanonical(Circle(),tocanonical(S,z)))
function evaluate{D<:Circle,R}(f::AbstractVector,S::Taylor{D,R},z)
    z=mappoint(S,𝕌,z)
    d=domain(S)
    horner(f,z)
end

function evaluate{D<:Domain}(f::AbstractVector,S::Hardy{false,D},z)
    z=mappoint(S,𝕌,z)
    z=1/z
    z*horner(f,z)
end
function evaluate{D<:Circle}(f::AbstractVector,S::Hardy{false,D},z)
    z=mappoint(S,𝕌,z)
    z=1/z
    z*horner(f,z)
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

points(sp::CosSpace,n) = points(domain(sp),2n-2)[1:n]




plan_transform(::CosSpace,x::AbstractVector) = plan_chebyshevtransform(x;kind=2)
plan_itransform(::CosSpace,x::AbstractVector) = plan_ichebyshevtransform(x;kind=2)
transform(::CosSpace,vals,plan) = plan*vals
itransform(::CosSpace,cfs,plan) = plan*cfs

evaluate(f::AbstractVector,S::CosSpace,t) = clenshaw(Chebyshev(),f,cos(tocanonical(S,t)))


points(sp::SinSpace,n)=points(domain(sp),2n+2)[2:n+1]

for (Typ,Pltr!,Pltr) in ((:TransformPlan,:plan_transform!,:plan_transform),
                         (:ITransformPlan,:plan_itransform!,:plan_itransform))
    @eval begin
        $Pltr!{DD,T<:FFTW.fftwNumber}(sp::SinSpace{DD},x::AbstractVector{T}) =
            $Typ(sp,FFTW.plan_r2r!(x,FFTW.RODFT00),Val{true})
        $Pltr{DD,T<:FFTW.fftwNumber}(sp::SinSpace{DD},x::AbstractVector{T}) =
            $Typ(sp,$Pltr!(sp,x),Val{false})
        $Pltr!{DD,T}(sp::SinSpace{DD},x::AbstractVector{T}) =
            error("transform for SinSpace only implemented for fftwNumbers")
        $Pltr{DD,T}(sp::SinSpace{DD},x::AbstractVector{T}) =
            error("transform for SinSpace only implemented for fftwNumbers")

        *{T,D,R}(P::$Typ{T,SinSpace{D,R},false},vals::AbstractVector{T}) = P.plan*copy(vals)
    end
end


*{T,DD,RR}(P::TransformPlan{T,SinSpace{DD,RR},true},vals::AbstractVector{T}) =
    scale!(one(T)/(length(vals)+1),P.plan*vals)
*{T,DD,RR}(P::ITransformPlan{T,SinSpace{DD,RR},true},cfs::AbstractVector{T}) =
    scale!(one(T)/2,P.plan*cfs)


transform(sp::SinSpace,vals::AbstractVector,plan) = plan*vals
itransform(sp::SinSpace,vals::AbstractVector,plan) = plan*vals

evaluate(f::AbstractVector,S::SinSpace,t) = sineshaw(f,tocanonical(S,t))



## Laurent space
doc"""
`Laurent()` is the space spanned by the complex exponentials
```
    1,exp(-im*θ),exp(im*θ),exp(-2im*θ),…
```
See also `Fourier`.
"""
const Laurent{DD,RR} = SumSpace{Tuple{Hardy{true,DD,RR},Hardy{false,DD,RR}},DD,RR}


##FFT That interlaces coefficients

plan_transform!{DD,RR,T<:Complex}(sp::Laurent{DD,RR},x::AbstractVector{T}) =
    TransformPlan(sp,plan_fft!(x),Val{true})
plan_itransform!{DD,RR,T<:Complex}(sp::Laurent{DD,RR},x::AbstractVector{T}) =
    ITransformPlan(sp,plan_ifft!(x),Val{true})

plan_transform!{DD,RR,T}(sp::Laurent{DD,RR},x::AbstractVector{T}) =
    error("In place variants not possible with real data.")
plan_itransform!{DD,RR,T}(sp::Laurent{DD,RR},x::AbstractVector{T}) =
    error("In place variants not possible with real data.")


plan_transform{T<:Complex,DD,RR}(sp::Laurent{DD,RR},x::AbstractVector{T}) =
    TransformPlan(sp,plan_transform!(sp,x),Val{false})
plan_itransform{T<:Complex,DD,RR}(sp::Laurent{DD,RR},x::AbstractVector{T}) =
    ITransformPlan(sp,plan_itransform!(sp,x),Val{false})

function plan_transform{T,DD,RR}(sp::Laurent{DD,RR},x::AbstractVector{T})
    plan = plan_transform(sp,Array{Complex{T}}(length(x))) # we can reuse vector in itransform
    TransformPlan{T,typeof(sp),false,typeof(plan)}(sp,plan)
end
function plan_itransform{T,DD,RR}(sp::Laurent{DD,RR},x::AbstractVector{T})
    plan = plan_itransform(sp,Array{Complex{T}}(length(x))) # we can reuse vector in itransform
    ITransformPlan{T,typeof(sp),false,typeof(plan)}(sp,plan)
end

function *{T,DD,RR}(P::TransformPlan{T,Laurent{DD,RR},true},vals::AbstractVector{T})
    n = length(vals)
    vals = scale!(inv(T(n)),P.plan*vals)
    reverseeven!(interlace!(vals,1))
end

function *{T,DD,RR}(P::ITransformPlan{T,Laurent{DD,RR},true},cfs::AbstractVector{T})
    n = length(cfs)
    reverseeven!(cfs)
    cfs[:]=[cfs[1:2:end];cfs[2:2:end]]  # TODO: deinterlace!
    scale!(n,cfs)
    P.plan*cfs
end

*{T<:Complex,DD,RR}(P::TransformPlan{T,Laurent{DD,RR},false},vals::AbstractVector{T}) = P.plan*copy(vals)
*{T,DD,RR}(P::TransformPlan{T,Laurent{DD,RR},false},vals::AbstractVector{T}) = P.plan*Vector{Complex{T}}(vals)
*{T<:Complex,DD,RR}(P::ITransformPlan{T,Laurent{DD,RR},false},vals::AbstractVector{T}) = P.plan*copy(vals)
*{T,DD,RR}(P::ITransformPlan{T,Laurent{DD,RR},false},vals::AbstractVector{T}) = P.plan*Vector{Complex{T}}(vals)



transform{DD,RR}(::Laurent{DD,RR},vals,plan) = plan*vals
itransform{DD,RR}(::Laurent{DD,RR},cfs,plan) = plan*cfs

transform{DD,RR}(sp::Laurent{DD,RR},vals::AbstractVector) = plan_transform(sp,vals)*vals
itransform{DD,RR}(sp::Laurent{DD,RR},cfs::AbstractVector) = plan_itransform(sp,cfs)*cfs






function evaluate{DD,RR}(f::AbstractVector,S::Laurent{DD,RR},z)
    z = mappoint(domain(S),Circle(),z)
    invz = 1./z
    horner(f,1:2:length(f),z) + horner(f,2:2:length(f),invz).*invz
end


function Base.conj{DD,RR}(f::Fun{Laurent{DD,RR}})
    cfs=Array{eltype(f)}(iseven(ncoefficients(f))?ncoefficients(f)+1:ncoefficients(f))
    cfs[1]=conj(f.coefficients[1])
    cfs[ncoefficients(f)] = 0
    for k=2:2:ncoefficients(f)-1
        cfs[k]=conj(f.coefficients[k+1])
    end
    for k=3:2:ncoefficients(f)+1
        cfs[k]=conj(f.coefficients[k-1])
    end
    Fun(space(f),cfs)
end

## Fourier space

doc"""
`Fourier()` is the space spanned by the trigonemtric polynomials
```
    1,sin(θ),cos(θ),sin(2θ),cos(2θ),…
```
See also `Laurent`.
"""
const Fourier{DD,RR} = SumSpace{Tuple{CosSpace{DD,RR},SinSpace{DD,RR}},DD,RR}

(::Type{Laurent})(d::Domain) = Laurent{typeof(d),complex(prectype(d))}(d)
(::Type{Fourier})(d::Domain) = Fourier{typeof(d),real(prectype(d))}(d)

for Typ in (:Laurent,:Fourier)
    @eval begin
        (::Type{$Typ})() = $Typ(PeriodicInterval())
        (::Type{$Typ})(d) = $Typ(PeriodicDomain(d))

        hasfasttransform{D,R}(::$Typ{D,R}) = true
    end
end


Laurent{DD,RR}(S::Fourier{DD,RR}) = Laurent(domain(S))
Fourier{DD,RR}(S::Laurent{DD,RR}) = Fourier(domain(S))

for T in (:CosSpace,:SinSpace)
    @eval begin
        # override default as canonicalspace must be implemented
        maxspace{D,R}(::$T,::Fourier{D,R}) = NoSpace()
        maxspace{D,R}(::Fourier{D,R},::$T) = NoSpace()
    end
end

points{D,R}(sp::Fourier{D,R},n)=points(domain(sp),n)

plan_transform!{T<:FFTW.fftwNumber,D,R}(sp::Fourier{D,R},x::AbstractVector{T}) =
    TransformPlan(sp,FFTW.plan_r2r!(x, FFTW.R2HC),Val{true})
plan_itransform!{T<:FFTW.fftwNumber,D,R}(sp::Fourier{D,R},x::AbstractVector{T}) =
    ITransformPlan(sp,FFTW.plan_r2r!(x, FFTW.HC2R),Val{true})

for (Typ,Pltr!,Pltr) in ((:TransformPlan,:plan_transform!,:plan_transform),
                         (:ITransformPlan,:plan_itransform!,:plan_itransform))
    @eval begin
        $Pltr{T<:FFTW.fftwNumber,DD,RR}(sp::Fourier{DD,RR},x::AbstractVector{T}) =
            $Typ(sp,$Pltr!(sp,x),Val{false})
        $Pltr!{T,DD,RR}(sp::Fourier{DD,RR},x::AbstractVector{T}) =
            error("transform for Fourier only implemented for fftwNumbers")
        $Pltr{T,DD,RR}(sp::Fourier{DD,RR},x::AbstractVector{T}) =
            error("transform for Fourier only implemented for fftwNumbers")

        *{T,DD,RR}(P::$Typ{T,Fourier{DD,RR},false},vals::AbstractVector{T}) = P.plan*copy(vals)
    end
end


function *{T,DD,RR}(P::TransformPlan{T,Fourier{DD,RR},true},vals::AbstractVector{T})
    n = length(vals)
    cfs = scale!(T(2)/n,P.plan*vals)
    cfs[1] /= 2
    if iseven(n)
        cfs[n÷2+1] /= 2
    end

    negateeven!(reverseeven!(interlace!(cfs,1)))
end

function *{T,DD,RR}(P::ITransformPlan{T,Fourier{DD,RR},true},cfs::AbstractVector{T})
    n = length(cfs)
    reverseeven!(negateeven!(cfs))
    cfs[:] = [cfs[1:2:end];cfs[2:2:end]]
    if iseven(n)
        cfs[n÷2+1] *= 2
    end
    cfs[1] *= 2
    P.plan*scale!(inv(T(2)),cfs)
end


transform{DD,RR}(sp::Fourier{DD,RR},vals::AbstractVector,plan) = plan*vals
itransform{DD,RR}(sp::Fourier{DD,RR},cfs::AbstractVector,plan) = plan*cfs

transform{DD,RR}(sp::Fourier{DD,RR},vals::AbstractVector) = plan_transform(sp,vals)*vals
itransform{DD,RR}(sp::Fourier{DD,RR},cfs::AbstractVector) = plan_itransform(sp,cfs)*cfs





canonicalspace{DD<:PeriodicInterval,RR}(S::Laurent{DD,RR}) = Fourier(domain(S))
canonicalspace{DD<:Circle,RR}(S::Fourier{DD,RR}) = Laurent(domain(S))
canonicalspace{DD<:PeriodicLine,RR}(S::Laurent{DD,RR}) = S


## Ones and zeros

for sp in (:Fourier,:CosSpace,:Laurent,:Taylor)
    @eval begin
        Base.ones{T<:Number,DD,RR}(::Type{T},S::$sp{DD,RR}) = Fun(S,ones(T,1))
        Base.ones{DD,RR}(S::$sp{DD,RR}) = Fun(S,ones(1))
    end
end


function identity_fun{DD<:Circle,RR}(S::Taylor{DD,RR})
    d=domain(S)
    if d.orientation
        Fun(S,[d.center,d.radius])
    else
        error("Cannot create identity on $S")
    end
end


identity_fun{DD<:Circle,RR}(S::Fourier{DD,RR}) = Fun(identity_fun(Laurent(domain(S))),S)


reverseorientation{DD,RR}(f::Fun{Fourier{DD,RR}}) =
    Fun(Fourier(reverse(domain(f))),alternatesign!(copy(f.coefficients)))
function reverseorientation{DD,RR}(f::Fun{Laurent{DD,RR}})
    # exp(im*k*x) -> exp(-im*k*x), or equivalentaly z -> 1/z
    n=ncoefficients(f)
    ret=Array{eltype(f)}(iseven(n)?n+1:n)  # since z -> 1/z we get one more coefficient
    ret[1]=f.coefficients[1]
    for k=2:2:length(ret)-1
        ret[k+1]=f.coefficients[k]
    end
    for k=2:2:n-1
        ret[k]=f.coefficients[k+1]
    end
    iseven(n) && (ret[n] = 0)

    Fun(Laurent(reverse(domain(f))),ret)
end

include("calculus.jl")
include("specialfunctions.jl")
include("FourierOperators.jl")
include("LaurentOperators.jl")
include("LaurentDirichlet.jl")
