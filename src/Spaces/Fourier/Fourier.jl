

export Fourier,Taylor,Hardy,CosSpace,SinSpace,Laurent

for T in (:CosSpace,:SinSpace)
    @eval begin
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
immutable Hardy{s,D<:Domain} <: UnivariateSpace{ComplexBasis,D}
    domain::D
    Hardy(d) = new(d)
    Hardy() = new(D())
end

# The <: Domain is crucial for matching Basecall overrides
typealias Taylor{D<:Domain} Hardy{true,D}


# Following is broken in 0.4
if VERSION ≥ v"0.5"
    doc"""
    `Taylor()` is the space spanned by `[1,z,z^2,...]`.  This is a type alias for `Hardy{true}`.

    """
    Taylor

    doc"""
    `Hardy{false}()` is the space spanned by `[1/z,1/z^2,...]`
    """
    Hardy{false}
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


for (Typ,Plfft!,Plfft,Pltr!,Pltr) in ((:TransformPlan,:plan_fft!,:plan_fft,:plan_transform!,:plan_transform),
                           (:ITransformPlan,:plan_ifft!,:plan_ifft,:plan_itransform!,:plan_itransform))
    @eval begin
        $Pltr!{T<:Complex}(sp::Hardy,x::Vector{T}) = $Typ(sp,$Plfft!(x),Val{true})
        $Pltr!{T<:Real}(::Hardy,x::Vector{T}) =
            error("In place variants not possible with real data.")

        $Pltr{T<:Complex}(sp::Hardy,x::Vector{T}) = $Typ(sp,$Pltr!(sp,x),Val{false})
        function $Pltr{T}(sp::Hardy,x::Vector{T})
            plan = $Pltr(sp,Array{Complex{T}}(length(x))) # we can reuse vector in itransform
            $Typ{T,typeof(sp),false,typeof(plan)}(sp,plan)
        end

        *{T<:Complex,HS<:Hardy}(P::$Typ{T,HS,false},vals::Vector{T}) = P.plan*copy(vals)
        *{T,HS<:Hardy}(P::$Typ{T,HS,false},vals::Vector{T}) = P.plan*Vector{Complex{T}}(vals)
    end
end


*{T,DD}(P::TransformPlan{T,Hardy{true,DD},true},vals::Vector{T}) =
    scale!(one(T)/length(vals),P.plan*vals)
*{T,DD}(P::ITransformPlan{T,Hardy{true,DD},true},cfs::Vector{T}) =
    scale!(length(cfs),P.plan*cfs)
*{T,DD}(P::TransformPlan{T,Hardy{false,DD},true},vals::Vector{T}) =
    scale!(one(T)/length(vals),reverse!(P.plan*vals))
*{T,DD}(P::ITransformPlan{T,Hardy{false,DD},true},cfs::Vector{T}) =
    scale!(length(cfs),P.plan*reverse!(cfs))


transform(sp::Hardy,vals::Vector,plan) = plan*vals
itransform(sp::Hardy,vals::Vector,plan) = plan*vals

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

points(sp::CosSpace,n) = points(domain(sp),2n-2)[1:n]




plan_transform(::CosSpace,x::Vector) = plan_chebyshevtransform(x;kind=2)
plan_itransform(::CosSpace,x::Vector) = plan_ichebyshevtransform(x;kind=2)
transform(::CosSpace,vals,plan) = plan*vals
itransform(::CosSpace,cfs,plan) = plan*cfs

evaluate(f::Vector,S::CosSpace,t) = clenshaw(Chebyshev(),f,cos(tocanonical(S,t)))


points(sp::SinSpace,n)=points(domain(sp),2n+2)[2:n+1]

for (Typ,Pltr!,Pltr) in ((:TransformPlan,:plan_transform!,:plan_transform),
                         (:ITransformPlan,:plan_itransform!,:plan_itransform))
    @eval begin
        $Pltr!{DD,T<:FFTW.fftwNumber}(sp::SinSpace{DD},x::Vector{T}) =
            $Typ(sp,FFTW.plan_r2r!(x,FFTW.RODFT00),Val{true})
        $Pltr{DD,T<:FFTW.fftwNumber}(sp::SinSpace{DD},x::Vector{T}) =
            $Typ(sp,$Pltr!(sp,x),Val{false})
        $Pltr!{DD,T}(sp::SinSpace{DD},x::Vector{T}) =
            error("transform for SinSpace only implemented for fftwNumbers")
        $Pltr{DD,T}(sp::SinSpace{DD},x::Vector{T}) =
            error("transform for SinSpace only implemented for fftwNumbers")

        *{T,D}(P::$Typ{T,SinSpace{D},false},vals::Vector{T}) = P.plan*copy(vals)
    end
end


*{T,DD}(P::TransformPlan{T,SinSpace{DD},true},vals::Vector{T}) =
    scale!(one(T)/(length(vals)+1),P.plan*vals)
*{T,DD}(P::ITransformPlan{T,SinSpace{DD},true},cfs::Vector{T}) =
    scale!(one(T)/2,P.plan*cfs)


transform(sp::SinSpace,vals::Vector,plan) = plan*vals
itransform(sp::SinSpace,vals::Vector,plan) = plan*vals

evaluate(f::AbstractVector,S::SinSpace,t) = sineshaw(f,tocanonical(S,t))



## Laurent space
doc"""
`Laurent()` is the space spanned by the complex exponentials
```
    1,exp(-im*θ),exp(im*θ),exp(-2im*θ),…
```
See also `Fourier`.
"""
typealias Laurent{DD} SumSpace{Tuple{Hardy{true,DD},Hardy{false,DD}},ComplexBasis,DD,1}


##FFT That interlaces coefficients

plan_transform!{DD,T<:Complex}(sp::Laurent{DD},x::Vector{T}) =
    TransformPlan(sp,plan_fft!(x),Val{true})
plan_itransform!{DD,T<:Complex}(sp::Laurent{DD},x::Vector{T}) =
    ITransformPlan(sp,plan_ifft!(x),Val{true})

plan_transform!{DD,T}(sp::Laurent{DD},x::Vector{T}) =
    error("In place variants not possible with real data.")
plan_itransform!{DD,T}(sp::Laurent{DD},x::Vector{T}) =
    error("In place variants not possible with real data.")


plan_transform{T<:Complex,DD}(sp::Laurent{DD},x::Vector{T}) =
    TransformPlan(sp,plan_transform!(sp,x),Val{false})
plan_itransform{T<:Complex,DD}(sp::Laurent{DD},x::Vector{T}) =
    ITransformPlan(sp,plan_itransform!(sp,x),Val{false})

function plan_transform{T,DD}(sp::Laurent{DD},x::Vector{T})
    plan = plan_transform(sp,Array{Complex{T}}(length(x))) # we can reuse vector in itransform
    TransformPlan{T,typeof(sp),false,typeof(plan)}(sp,plan)
end
function plan_itransform{T,DD}(sp::Laurent{DD},x::Vector{T})
    plan = plan_itransform(sp,Array{Complex{T}}(length(x))) # we can reuse vector in itransform
    ITransformPlan{T,typeof(sp),false,typeof(plan)}(sp,plan)
end

function *{T,DD}(P::TransformPlan{T,Laurent{DD},true},vals::Vector{T})
    n = length(vals)
    vals = scale!(inv(T(n)),P.plan*vals)
    reverseeven!(interlace!(vals,1))
end

function *{T,DD}(P::ITransformPlan{T,Laurent{DD},true},cfs::Vector{T})
    n = length(cfs)
    reverseeven!(cfs)
    cfs[:]=[cfs[1:2:end];cfs[2:2:end]]  # TODO: deinterlace!
    scale!(n,cfs)
    P.plan*cfs
end

*{T<:Complex,DD}(P::TransformPlan{T,Laurent{DD},false},vals::Vector{T}) = P.plan*copy(vals)
*{T,DD}(P::TransformPlan{T,Laurent{DD},false},vals::Vector{T}) = P.plan*Vector{Complex{T}}(vals)
*{T<:Complex,DD}(P::ITransformPlan{T,Laurent{DD},false},vals::Vector{T}) = P.plan*copy(vals)
*{T,DD}(P::ITransformPlan{T,Laurent{DD},false},vals::Vector{T}) = P.plan*Vector{Complex{T}}(vals)



transform{DD}(::Laurent{DD},vals,plan) = plan*vals
itransform{DD}(::Laurent{DD},cfs,plan) = plan*cfs

transform{DD}(sp::Laurent{DD},vals::Vector) = plan_transform(sp,vals)*vals
itransform{DD}(sp::Laurent{DD},cfs::Vector) = plan_itransform(sp,cfs)*cfs






function evaluate{DD}(f::AbstractVector,S::Laurent{DD},z)
    z = mappoint(domain(S),Circle(),z)
    invz = 1./z
    horner(f,1:2:length(f),z) + horner(f,2:2:length(f),invz).*invz
end


function Base.conj{DD}(f::Fun{Laurent{DD}})
    cfs=Array(eltype(f),iseven(ncoefficients(f))?ncoefficients(f)+1:ncoefficients(f))
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
typealias Fourier{DD} SumSpace{Tuple{CosSpace{DD},SinSpace{DD}},RealBasis,DD,1}

for Typ in (:Laurent,:Fourier)
    @eval begin
        @compat (::Type{$Typ})(d::Domain) = $Typ{typeof(d)}(d)
        @compat (::Type{$Typ})() = $Typ(PeriodicInterval())
        @compat (::Type{$Typ})(d) = $Typ(PeriodicDomain(d))

        hasfasttransform{D}(::$Typ{D}) = true
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

plan_transform!{T<:FFTW.fftwNumber,D}(sp::Fourier{D},x::Vector{T}) =
    TransformPlan(sp,FFTW.plan_r2r!(x, FFTW.R2HC),Val{true})
plan_itransform!{T<:FFTW.fftwNumber,D}(sp::Fourier{D},x::Vector{T}) =
    TransformPlan(sp,FFTW.plan_r2r!(x, FFTW.HC2R),Val{true})

for (Typ,Pltr!,Pltr) in ((:TransformPlan,:plan_transform!,:plan_transform),
                         (:ITransformPlan,:plan_itransform!,:plan_itransform))
    @eval begin
        $Pltr{T<:FFTW.fftwNumber,DD}(sp::Fourier{DD},x::Vector{T}) =
            $Typ(sp,$Pltr!(sp,x),Val{false})
        $Pltr!{T,DD}(sp::Fourier{DD},x::Vector{T}) =
            error("transform for Fourier only implemented for fftwNumbers")
        $Pltr{T,DD}(sp::Fourier{DD},x::Vector{T}) =
            error("transform for Fourier only implemented for fftwNumbers")

        *{T,DD}(P::$Typ{T,Fourier{DD},false},vals::Vector{T}) = P.plan*copy(vals)
    end
end


function *{T,DD}(P::TransformPlan{T,Fourier{DD},true},vals::Vector{T})
    n = length(vals)
    cfs = scale!(T(2)/n,P.plan*vals)
    cfs[1] /= 2
    if iseven(n)
        cfs[n÷2+1] /= 2
    end

    negateeven!(reverseeven!(interlace!(cfs,1)))
end

function *{T,DD}(P::ITransformPlan{T,Fourier{DD},true},cfs::Vector{T})
    n = length(cfs)
    reverseeven!(negativeeven!(cfs))
    cfs[:] = [cfs[1:2:end];cfs[2:2:end]]
    if iseven(n)
        cfs[n÷2+1] *= 2
    end
    cfs[1] *= 2
    P.plan*scale!(inv(T(2)),cfs)
end


transform{DD}(sp::Fourier{DD},vals::Vector,plan) = plan*vals
itransform{DD}(sp::Fourier{DD},cfs::Vector,plan) = plan*cfs

transform{DD}(sp::Fourier{DD},vals::Vector) = plan_transform(sp,vals)*vals
itransform{DD}(sp::Fourier{DD},cfs::Vector) = plan_itransform(sp,cfs)*cfs





canonicalspace{DD<:PeriodicInterval}(S::Laurent{DD})=Fourier(domain(S))
canonicalspace{DD<:Circle}(S::Fourier{DD})=Laurent(domain(S))
canonicalspace{DD<:PeriodicLine}(S::Laurent{DD})=S

for Typ in (:CosSpace,:Taylor)
    @eval union_rule(A::ConstantSpace,B::$Typ)=B
end

## Ones and zeros

for sp in (:Fourier,:CosSpace,:Laurent,:Taylor)
    @eval begin
        Base.ones{T<:Number,D}(::Type{T},S::$sp{D})=Fun(S,ones(T,1))
        Base.ones{D}(S::$sp{D})=Fun(S,ones(1))
    end
end


function identity_fun{DD<:Circle}(S::Taylor{DD})
    d=domain(S)
    if d.orientation
        Fun(S,[d.center,d.radius])
    else
        error("Cannot create identity on $S")
    end
end


identity_fun{DD<:Circle}(S::Fourier{DD}) = Fun(identity_fun(Laurent(domain(S))),S)


reverseorientation{D}(f::Fun{Fourier{D}}) =
    Fun(Fourier(reverse(domain(f))),alternatesign!(copy(f.coefficients)))
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

    Fun(Laurent(reverse(domain(f))),ret)
end

include("calculus.jl")
include("specialfunctions.jl")
include("FourierOperators.jl")
include("LaurentOperators.jl")
include("LaurentDirichlet.jl")
