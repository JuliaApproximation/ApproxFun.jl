

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



Base.promote_rule(::Type{Fun{S,V}},::Type{T}) where {T<:Number,S<:Union{Hardy{true},CosSpace},V} =
    Fun{S,promote_type(V,T)}
Base.promote_rule(::Type{Fun{S}},::Type{T}) where {T<:Number,S<:Union{Hardy{true},CosSpace}} =
    Fun{S,T}

(H::Type{Hardy{s}})(d::Domain) where {s} = Hardy{s,typeof(d),complex(prectype(d))}(d)
(H::Type{Hardy{s}})() where {s} = Hardy{s}(Circle())

canonicalspace(S::Hardy) = S
setdomain(S::Hardy{s},d::Domain) where {s} = Hardy{s}(d)


spacescompatible(a::Hardy{s},b::Hardy{s}) where {s} = domainscompatible(a,b)
hasfasttransform(::Hardy) = true


for (Typ,Plfft!,Plfft,Pltr!,Pltr) in ((:TransformPlan,:plan_fft!,:plan_fft,:plan_transform!,:plan_transform),
                           (:ITransformPlan,:plan_ifft!,:plan_ifft,:plan_itransform!,:plan_itransform))
    @eval begin
        $Pltr!(sp::Hardy,x::AbstractVector{T}) where {T<:Complex} = $Typ(sp,$Plfft!(x),Val{true})
        $Pltr!(::Hardy,x::AbstractVector{T}) where {T<:Real} =
            error("In place variants not possible with real data.")

        $Pltr(sp::Hardy,x::AbstractVector{T}) where {T<:Complex} = $Typ(sp,$Pltr!(sp,x),Val{false})
        function $Pltr(sp::Hardy,x::AbstractVector{T}) where T
            plan = $Pltr(sp,Array{Complex{T}}(length(x))) # we can reuse vector in itransform
            $Typ{T,typeof(sp),false,typeof(plan)}(sp,plan)
        end

        *(P::$Typ{T,HS,false},vals::AbstractVector{T}) where {T<:Complex,HS<:Hardy} = P.plan*copy(vals)
        *(P::$Typ{T,HS,false},vals::AbstractVector{T}) where {T,HS<:Hardy} = P.plan*Vector{Complex{T}}(vals)
    end
end


*(P::TransformPlan{T,Hardy{true,DD,RR},true},vals::AbstractVector{T}) where {T,DD,RR} =
    scale!(one(T)/length(vals),P.plan*vals)
*(P::ITransformPlan{T,Hardy{true,DD,RR},true},cfs::AbstractVector{T}) where {T,DD,RR} =
    scale!(length(cfs),P.plan*cfs)
*(P::TransformPlan{T,Hardy{false,DD,RR},true},vals::AbstractVector{T}) where {T,DD,RR} =
    scale!(one(T)/length(vals),reverse!(P.plan*vals))
*(P::ITransformPlan{T,Hardy{false,DD,RR},true},cfs::AbstractVector{T}) where {T,DD,RR} =
    scale!(length(cfs),P.plan*reverse!(cfs))


transform(sp::Hardy,vals::AbstractVector,plan) = plan*vals
itransform(sp::Hardy,vals::AbstractVector,plan) = plan*vals

evaluate(f::AbstractVector,S::Taylor{D,R},z) where {D<:Domain,R} = horner(f,fromcanonical(Circle(),tocanonical(S,z)))
function evaluate(f::AbstractVector,S::Taylor{D,R},z) where {D<:Circle,R}
    z=mappoint(S,𝕌,z)
    d=domain(S)
    horner(f,z)
end

function evaluate(f::AbstractVector,S::Hardy{false,D},z) where D<:Domain
    z=mappoint(S,𝕌,z)
    z=1/z
    z*horner(f,z)
end
function evaluate(f::AbstractVector,S::Hardy{false,D},z) where D<:Circle
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
        $Pltr!(sp::SinSpace{DD},x::AbstractVector{T}) where {DD,T<:FFTW.fftwNumber} =
            $Typ(sp,FFTW.plan_r2r!(x,FFTW.RODFT00),Val{true})
        $Pltr(sp::SinSpace{DD},x::AbstractVector{T}) where {DD,T<:FFTW.fftwNumber} =
            $Typ(sp,$Pltr!(sp,x),Val{false})
        $Pltr!(sp::SinSpace{DD},x::AbstractVector{T}) where {DD,T} =
            error("transform for SinSpace only implemented for fftwNumbers")
        $Pltr(sp::SinSpace{DD},x::AbstractVector{T}) where {DD,T} =
            error("transform for SinSpace only implemented for fftwNumbers")

        *(P::$Typ{T,SinSpace{D,R},false},vals::AbstractVector{T}) where {T,D,R} = P.plan*copy(vals)
    end
end


*(P::TransformPlan{T,SinSpace{DD,RR},true},vals::AbstractVector{T}) where {T,DD,RR} =
    scale!(one(T)/(length(vals)+1),P.plan*vals)
*(P::ITransformPlan{T,SinSpace{DD,RR},true},cfs::AbstractVector{T}) where {T,DD,RR} =
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

plan_transform!(sp::Laurent{DD,RR},x::AbstractVector{T}) where {DD,RR,T<:Complex} =
    TransformPlan(sp,plan_fft!(x),Val{true})
plan_itransform!(sp::Laurent{DD,RR},x::AbstractVector{T}) where {DD,RR,T<:Complex} =
    ITransformPlan(sp,plan_ifft!(x),Val{true})

plan_transform!(sp::Laurent{DD,RR},x::AbstractVector{T}) where {DD,RR,T} =
    error("In place variants not possible with real data.")
plan_itransform!(sp::Laurent{DD,RR},x::AbstractVector{T}) where {DD,RR,T} =
    error("In place variants not possible with real data.")


plan_transform(sp::Laurent{DD,RR},x::AbstractVector{T}) where {T<:Complex,DD,RR} =
    TransformPlan(sp,plan_transform!(sp,x),Val{false})
plan_itransform(sp::Laurent{DD,RR},x::AbstractVector{T}) where {T<:Complex,DD,RR} =
    ITransformPlan(sp,plan_itransform!(sp,x),Val{false})

function plan_transform(sp::Laurent{DD,RR},x::AbstractVector{T}) where {T,DD,RR}
    plan = plan_transform(sp,Array{Complex{T}}(length(x))) # we can reuse vector in itransform
    TransformPlan{T,typeof(sp),false,typeof(plan)}(sp,plan)
end
function plan_itransform(sp::Laurent{DD,RR},x::AbstractVector{T}) where {T,DD,RR}
    plan = plan_itransform(sp,Array{Complex{T}}(length(x))) # we can reuse vector in itransform
    ITransformPlan{T,typeof(sp),false,typeof(plan)}(sp,plan)
end

function *(P::TransformPlan{T,Laurent{DD,RR},true},vals::AbstractVector{T}) where {T,DD,RR}
    n = length(vals)
    vals = scale!(inv(T(n)),P.plan*vals)
    reverseeven!(interlace!(vals,1))
end

function *(P::ITransformPlan{T,Laurent{DD,RR},true},cfs::AbstractVector{T}) where {T,DD,RR}
    n = length(cfs)
    reverseeven!(cfs)
    cfs[:]=[cfs[1:2:end];cfs[2:2:end]]  # TODO: deinterlace!
    scale!(n,cfs)
    P.plan*cfs
end

*(P::TransformPlan{T,Laurent{DD,RR},false},vals::AbstractVector{T}) where {T<:Complex,DD,RR} = P.plan*copy(vals)
*(P::TransformPlan{T,Laurent{DD,RR},false},vals::AbstractVector{T}) where {T,DD,RR} = P.plan*Vector{Complex{T}}(vals)
*(P::ITransformPlan{T,Laurent{DD,RR},false},vals::AbstractVector{T}) where {T<:Complex,DD,RR} = P.plan*copy(vals)
*(P::ITransformPlan{T,Laurent{DD,RR},false},vals::AbstractVector{T}) where {T,DD,RR} = P.plan*Vector{Complex{T}}(vals)



transform(::Laurent{DD,RR},vals,plan) where {DD,RR} = plan*vals
itransform(::Laurent{DD,RR},cfs,plan) where {DD,RR} = plan*cfs

transform(sp::Laurent{DD,RR},vals::AbstractVector) where {DD,RR} = plan_transform(sp,vals)*vals
itransform(sp::Laurent{DD,RR},cfs::AbstractVector) where {DD,RR} = plan_itransform(sp,cfs)*cfs






function evaluate(f::AbstractVector,S::Laurent{DD,RR},z) where {DD,RR}
    z = mappoint(domain(S),Circle(),z)
    invz = 1./z
    horner(f,1:2:length(f),z) + horner(f,2:2:length(f),invz).*invz
end


function Base.conj(f::Fun{Laurent{DD,RR}}) where {DD,RR}
    ncoefficients(f) == 0 && return f

    cfs = Array{eltype(f)}(iseven(ncoefficients(f))?ncoefficients(f)+1:ncoefficients(f))
    cfs[1] = conj(f.coefficients[1])
    cfs[ncoefficients(f)] = 0
    for k=2:2:ncoefficients(f)-1
        cfs[k] = conj(f.coefficients[k+1])
    end
    for k=3:2:ncoefficients(f)+1
        cfs[k] = conj(f.coefficients[k-1])
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

Laurent(d::Domain) = Laurent{typeof(d),complex(prectype(d))}(d)
Fourier(d::Domain) = Fourier{typeof(d),real(prectype(d))}(d)

for Typ in (:Laurent,:Fourier)
    @eval begin
        (::Type{$Typ})() = $Typ(PeriodicInterval())
        (::Type{$Typ})(d) = $Typ(PeriodicDomain(d))

        hasfasttransform(::$Typ{D,R}) where {D,R} = true
    end
end


Laurent(S::Fourier{DD,RR}) where {DD,RR} = Laurent(domain(S))
Fourier(S::Laurent{DD,RR}) where {DD,RR} = Fourier(domain(S))

for T in (:CosSpace,:SinSpace)
    @eval begin
        # override default as canonicalspace must be implemented
        maxspace(::$T,::Fourier{D,R}) where {D,R} = NoSpace()
        maxspace(::Fourier{D,R},::$T) where {D,R} = NoSpace()
    end
end

points(sp::Fourier{D,R},n) where {D,R}=points(domain(sp),n)

plan_transform!(sp::Fourier{D,R},x::AbstractVector{T}) where {T<:FFTW.fftwNumber,D,R} =
    TransformPlan(sp,FFTW.plan_r2r!(x, FFTW.R2HC),Val{true})
plan_itransform!(sp::Fourier{D,R},x::AbstractVector{T}) where {T<:FFTW.fftwNumber,D,R} =
    ITransformPlan(sp,FFTW.plan_r2r!(x, FFTW.HC2R),Val{true})

for (Typ,Pltr!,Pltr) in ((:TransformPlan,:plan_transform!,:plan_transform),
                         (:ITransformPlan,:plan_itransform!,:plan_itransform))
    @eval begin
        $Pltr(sp::Fourier{DD,RR},x::AbstractVector{T}) where {T<:FFTW.fftwNumber,DD,RR} =
            $Typ(sp,$Pltr!(sp,x),Val{false})
        $Pltr!(sp::Fourier{DD,RR},x::AbstractVector{T}) where {T,DD,RR} =
            error("transform for Fourier only implemented for fftwNumbers")
        $Pltr(sp::Fourier{DD,RR},x::AbstractVector{T}) where {T,DD,RR} =
            error("transform for Fourier only implemented for fftwNumbers")

        *(P::$Typ{T,Fourier{DD,RR},false},vals::AbstractVector{T}) where {T,DD,RR} = P.plan*copy(vals)
    end
end


function *(P::TransformPlan{T,Fourier{DD,RR},true},vals::AbstractVector{T}) where {T,DD,RR}
    n = length(vals)
    cfs = scale!(T(2)/n,P.plan*vals)
    cfs[1] /= 2
    if iseven(n)
        cfs[n÷2+1] /= 2
    end

    negateeven!(reverseeven!(interlace!(cfs,1)))
end

function *(P::ITransformPlan{T,Fourier{DD,RR},true},cfs::AbstractVector{T}) where {T,DD,RR}
    n = length(cfs)
    reverseeven!(negateeven!(cfs))
    cfs[:] = [cfs[1:2:end];cfs[2:2:end]]
    if iseven(n)
        cfs[n÷2+1] *= 2
    end
    cfs[1] *= 2
    P.plan*scale!(inv(T(2)),cfs)
end


transform(sp::Fourier{DD,RR},vals::AbstractVector,plan) where {DD,RR} = plan*vals
itransform(sp::Fourier{DD,RR},cfs::AbstractVector,plan) where {DD,RR} = plan*cfs

transform(sp::Fourier{DD,RR},vals::AbstractVector) where {DD,RR} = plan_transform(sp,vals)*vals
itransform(sp::Fourier{DD,RR},cfs::AbstractVector) where {DD,RR} = plan_itransform(sp,cfs)*cfs





canonicalspace(S::Laurent{DD,RR}) where {DD<:PeriodicInterval,RR} = Fourier(domain(S))
canonicalspace(S::Fourier{DD,RR}) where {DD<:Circle,RR} = Laurent(domain(S))
canonicalspace(S::Laurent{DD,RR}) where {DD<:PeriodicLine,RR} = S


## Ones and zeros

for sp in (:Fourier,:CosSpace,:Laurent,:Taylor)
    @eval begin
        Base.ones(::Type{T},S::$sp{DD,RR}) where {T<:Number,DD,RR} = Fun(S,ones(T,1))
        Base.ones(S::$sp{DD,RR}) where {DD,RR} = Fun(S,ones(1))
    end
end


function identity_fun(S::Taylor{DD,RR}) where {DD<:Circle,RR}
    d=domain(S)
    if d.orientation
        Fun(S,[d.center,d.radius])
    else
        error("Cannot create identity on $S")
    end
end


identity_fun(S::Fourier{DD,RR}) where {DD<:Circle,RR} = Fun(identity_fun(Laurent(domain(S))),S)


reverseorientation(f::Fun{Fourier{DD,RR}}) where {DD,RR} =
    Fun(Fourier(reverse(domain(f))),alternatesign!(copy(f.coefficients)))
function reverseorientation(f::Fun{Laurent{DD,RR}}) where {DD,RR}
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
