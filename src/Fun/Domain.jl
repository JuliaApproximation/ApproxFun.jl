

export Domain,IntervalDomain,PeriodicDomain,tocanonical,fromcanonical,fromcanonicalD,∂
export chebyshevpoints,fourierpoints,isambiguous,arclength


# T is the numeric type used to represent the domain
# For d-dimensional domains, it is Vec{d,T}

abstract type Domain{T} end
const UnivariateDomain{T} = Domain{T} where {T<:Number}
const BivariateDomain{T} = Domain{Vec{2,T}} where {T<:Number}


eltype{T}(::Domain{T}) = T
eltype{T}(::Type{Domain{T}}) = T
Base.isreal{T<:Real}(::Domain{T}) = true
Base.isreal{T}(::Domain{T}) = false
dimension(::Domain{<:Number}) = 1
dimension(::Type{Domain{Vec{d,T}}}) where {T,d} = d
dimension(::Type{DT}) where {DT<:Domain} = dimension(supertype(DT))

# add indexing for all spaces, not just DirectSumSpace
# mimicking scalar vs vector

# TODO: 0.5 iterator
Base.start(s::Domain) = false
Base.next(s::Domain,st) = (s,true)
Base.done(s::Domain,st) = st
Base.length(s::Domain) = 1
getindex(s::Domain,::CartesianIndex{0}) = s
getindex(s::Domain,k) = k == 1 ? s : throw(BoundsError())
Base.endof(s::Domain) = 1


#supports broadcasting, overloaded for ArraySpace
Base.size(::Domain) = ()


# prectype gives the precision, including for Vec
prectype(d::Domain) = eltype(eltype(d))
prectype(::Type{D}) where {D<:Domain} = eltype(eltype(D))

#TODO: bivariate AnyDomain
struct AnyDomain <: Domain{UnsetNumber} end
struct EmptyDomain <: Domain{Void} end

isambiguous(::AnyDomain) = true
dimension(::AnyDomain) = 1

complexlength(::AnyDomain) = NaN
arclength(::AnyDomain) = NaN

Base.reverse(a::Union{AnyDomain,EmptyDomain}) = a

canonicaldomain(a::Union{AnyDomain,EmptyDomain}) = a

Base.in(x::Domain,::EmptyDomain) = false

##General routines


Base.isempty(::EmptyDomain) = true
Base.isempty(::Domain) = false
Base.intersect(a::Domain,b::Domain) = a==b ? a : EmptyDomain()


# TODO: throw error for override
Base.setdiff(a::Domain,b) = a == b ? EmptyDomain() : a
\(a::Domain,b) = setdiff(a,b)

## Interval Domains

abstract type IntervalDomain{T} <: UnivariateDomain{T} end

canonicaldomain(d::IntervalDomain) = Segment{real(prectype(d))}()

Base.isapprox(a::Domain,b::Domain) = a==b
domainscompatible(a,b) = domainscompatible(domain(a),domain(b))
domainscompatible(a::Domain,b::Domain) = isambiguous(a) || isambiguous(b) ||
                    isapprox(a,b)

function chebyshevpoints{T<:Number}(::Type{T},n::Integer;kind::Integer=1)
    if kind == 1
        T[sinpi((n-2k-one(T))/2n) for k=0:n-1]
    elseif kind == 2
        if n==1
            zeros(T,1)
        else
            T[cospi(k/(n-one(T))) for k=0:n-1]
        end
    end
end
chebyshevpoints(n::Integer;kind::Integer=1) = chebyshevpoints(Float64,n;kind=kind)

##TODO: Should fromcanonical be fromcanonical!?

points{T}(d::IntervalDomain{T},n::Integer) =
    fromcanonical.(d,chebyshevpoints(real(eltype(T)),n))  # eltype to handle point
bary(v::AbstractVector{Float64},d::IntervalDomain,x::Float64) = bary(v,tocanonical(d,x))

#TODO consider moving these
Base.first{T}(d::IntervalDomain{T}) = fromcanonical(d,-one(T))
Base.last{T}(d::IntervalDomain{T}) = fromcanonical(d,one(T))

Base.in(x,::AnyDomain) = true
function Base.in(x,d::IntervalDomain)
    T=real(prectype(d))
    y=tocanonical(d,x)
    ry=real(y)
    iy=imag(y)
    sc=norm(fromcanonicalD(d,ry<-1?-1:(ry>1?1:ry)))  # scale based on stretch of map on projection to interal
    isapprox(fromcanonical(d,y),x) &&
        -one(T)-100eps(T)/sc≤ry≤one(T)+100eps(T)/sc &&
        -100eps(T)/sc≤iy≤100eps(T)/sc
end

components(d::Domain) = [d]
issubcomponent(a::Domain,b::Domain) = a in components(b)

###### Periodic domains

abstract type PeriodicDomain{T} <: UnivariateDomain{T} end


canonicaldomain(::PeriodicDomain)=PeriodicInterval()

points{T}(d::PeriodicDomain{T},n::Integer) =
    fromcanonical.(d, fourierpoints(real(eltype(T)),n))

fourierpoints(n::Integer) = fourierpoints(Float64,n)
fourierpoints{T<:Number}(::Type{T},n::Integer)= convert(T,π)*collect(0:2:2n-2)/n


function Base.in{T}(x,d::PeriodicDomain{T})
    y=tocanonical(d,x)
    if !isapprox(fromcanonical(d,y),x)
        return false
    end

    l=arclength(d)
    if isinf(l)
        abs(imag(y))<20eps(T) && -2eps(T)<real(y)<2π+2eps(T)
    else
        abs(imag(y))/l<20eps(T) && -2l*eps(T)<real(y)<2π+2l*eps(T)
    end
end

Base.issubset(a::Domain,b::Domain)=a==b


Base.first(d::PeriodicDomain) = fromcanonical(d,0)
Base.last(d::PeriodicDomain) = fromcanonical(d,2π)


struct AnyPeriodicDomain <: PeriodicDomain{UnsetNumber} end
isambiguous(::AnyPeriodicDomain)=true

Base.convert{D<:PeriodicDomain}(::Type{D},::AnyDomain)=AnyPeriodicDomain()

## conveninece routines

Base.ones(d::Domain)=ones(prectype(d),Space(d))
Base.zeros(d::Domain)=zeros(prectype(d),Space(d))





function commondomain(P::AbstractVector)
    ret = AnyDomain()

    for op in P
        d = domain(op)
        @assert ret == AnyDomain() || d == AnyDomain() || ret == d

        if d != AnyDomain()
            ret = d
        end
    end

    ret
end

commondomain{T<:Number}(P::AbstractVector,g::Array{T}) = commondomain(P)
commondomain(P::AbstractVector,g) = commondomain([P;g])


domain(::Number) = AnyDomain()




## rand


Base.rand(d::IntervalDomain,k...) = fromcanonical.(d,2rand(k...)-1)
Base.rand(d::PeriodicDomain,k...) = fromcanonical.(d,2π*rand(k...)-π)

checkpoints(d::IntervalDomain) = fromcanonical.(d,[-0.823972,0.01,0.3273484])
checkpoints(d::PeriodicDomain) = fromcanonical.(d,[1.223972,3.14,5.83273484])

## boundary

doc"""
    ∂(d::Domain)

returns the boundary of `d`.  For example, the boundary of a `Disk()`
is a `Circle()`, and the boundary of `Interval()^2` is a piecewise interval
that sketches the boundary of a rectangle.
"""
∂(d::Domain) = EmptyDomain()   # This is meant to be overriden
∂(d::IntervalDomain) = [first(d),last(d)] #TODO: Points domain
∂(d::PeriodicDomain) = EmptyDomain()




## map domains
# we auto vectorize arguments
tocanonical(d::Domain,x,y,z...) = tocanonical(d,Vec(x,y,z...))
fromcanonical(d::Domain,x,y,z...) = fromcanonical(d,Vec(x,y,z...))


mappoint(d1::Domain,d2::Domain,x...) = fromcanonical(d2,tocanonical(d1,x...))
invfromcanonicalD(d::Domain,x...) = 1./fromcanonicalD(d,x...)



## domains in higher dimensions


## sorting
# we sort spaces lexigraphically by default

for OP in (:<,:(<=),:>,:(>=),:(Base.isless))
    @eval $OP(a::Domain,b::Domain)=$OP(string(a),string(b))
end


## Other special domains

struct PositiveIntegers <: Domain{Int} end
struct Integers <: Domain{Int} end

const ℕ = PositiveIntegers()
const ℤ = Integers()
