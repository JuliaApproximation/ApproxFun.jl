

export Domain, IntervalDomain, PeriodicDomain, tocanonical, fromcanonical, fromcanonicalD, ∂
export chebyshevpoints, fourierpoints, isambiguous, arclength
export components, component, ncomponents


# T is the numeric type used to represent the domain
# For d-dimensional domains, it is Vec{d,T}

abstract type Domain{T} end
const UnivariateDomain{T} = Domain{T} where {T<:Number}
const BivariateDomain{T} = Domain{Vec{2,T}} where {T<:Number}

struct DomainStyle <: BroadcastStyle end

BroadcastStyle(::Type{<:Domain}) = DomainStyle()
BroadcastStyle(A::AbstractArrayStyle, ::DomainStyle)  = A
BroadcastStyle(::DomainStyle, A::AbstractArrayStyle)  = A

Base.broadcast_axes(::Type{<:Domain}, A) = axes(A)
Base.broadcastable(x::Domain) = x

eltype(::Domain{T}) where T = T
eltype(::Type{D}) where D<:Domain{T} where T = T
isreal(::Domain{T}) where T<:Real = true
isreal(::Domain{T}) where T = false

copy(d::Domain) = d  # all domains are immutable

dimension(::Type{Domain{TT}}) where TT<:Number = 1
dimension(::Type{Domain{Vec{d,T}}}) where {T,d} = d
dimension(::Type{DT}) where {DT<:Domain} = dimension(supertype(DT))

dimension(d::Domain) = dimension(typeof(d))

# add indexing for all spaces, not just DirectSumSpace
# mimicking scalar vs vector

# prectype gives the precision, including for Vec
prectype(d::Domain) = eltype(eltype(d))
prectype(::Type{D}) where {D<:Domain} = eltype(eltype(D))

#TODO: bivariate AnyDomain
struct AnyDomain <: Domain{UnsetNumber} end
struct EmptyDomain <: Domain{Nothing} end

isambiguous(::AnyDomain) = true
dimension(::AnyDomain) = 1

complexlength(::AnyDomain) = NaN
arclength(::AnyDomain) = NaN

reverse(a::Union{AnyDomain,EmptyDomain}) = a

canonicaldomain(a::Union{AnyDomain,EmptyDomain}) = a

in(x::Domain,::EmptyDomain) = false

##General routines


isempty(::EmptyDomain) = true
isempty(::Domain) = false
intersect(a::Domain,b::Domain) = a==b ? a : EmptyDomain()


# TODO: throw error for override
setdiff(a::Domain,b) = a == b ? EmptyDomain() : a
\(a::Domain,b) = setdiff(a,b)

## Interval Domains

abstract type IntervalDomain{T} <: UnivariateDomain{T} end

canonicaldomain(d::IntervalDomain) = Segment{real(prectype(d))}()

isapprox(a::Domain,b::Domain) = a==b
domainscompatible(a,b) = domainscompatible(domain(a),domain(b))
domainscompatible(a::Domain,b::Domain) = isambiguous(a) || isambiguous(b) ||
                    isapprox(a,b)

##TODO: Should fromcanonical be fromcanonical!?

points(d::IntervalDomain{T},n::Integer;kind::Int=1) where {T} =
    fromcanonical.(Ref(d), chebyshevpoints(real(eltype(T)),n;kind=kind))  # eltype to handle point
bary(v::AbstractVector{Float64},d::IntervalDomain,x::Float64) = bary(v,tocanonical(d,x))

#TODO consider moving these
first(d::IntervalDomain{T}) where {T} = fromcanonical(d,-one(T))
last(d::IntervalDomain{T}) where {T} = fromcanonical(d,one(T))

in(x,::AnyDomain) = true
function in(x,d::IntervalDomain)
    T=float(real(prectype(d)))
    y=tocanonical(d,x)
    ry=real(y)
    iy=imag(y)
    sc=norm(fromcanonicalD(d,ry<-1 ? -one(ry) : (ry>1 ? one(ry) : ry)))  # scale based on stretch of map on projection to interal
    dy=fromcanonical(d,y)
    # TODO: use isapprox once keywords are fast
    ((isinf(norm(dy)) && isinf(norm(x))) ||  norm(dy-x) ≤ 1000eps(T)*max(norm(x),1)) &&
        -one(T)-100eps(T)/sc ≤ ry ≤ one(T)+100eps(T)/sc &&
        -100eps(T)/sc ≤ iy ≤ 100eps(T)/sc
end

ncomponents(s::Domain) = 1
components(s::Domain) = [s]
function components(s::Domain,k)
    k ≠ 1 && throw(BoundsError())
    s
end

issubcomponent(a::Domain,b::Domain) = a in components(b)


##### canoncial
"""
    canonicaldomain(d)

returns a domain which we map to for operations. For example,
the canonical domain for an interval [a,b] is [-1,1]
"""
function canonicaldomain end


"""
    tocanonical(d, x)

maps the point `x` in `d` to a point in `canonical(d,x)`
"""
function tocanonical end

###### Periodic domains

abstract type PeriodicDomain{T} <: UnivariateDomain{T} end


canonicaldomain(::PeriodicDomain) = PeriodicInterval()


points(d::PeriodicDomain{T},n::Integer) where {T} =
    fromcanonical.(Ref(d), fourierpoints(real(eltype(T)),n))

fourierpoints(n::Integer) = fourierpoints(Float64,n)
fourierpoints(::Type{T},n::Integer) where {T<:Number} = convert(T,π)*collect(0:2:2n-2)/n


function in(x,d::PeriodicDomain{T}) where T
    y=tocanonical(d,x)
    if !isapprox(fromcanonical(d,y),x)
        return false
    end

    l=arclength(d)
    if isinf(l)
        abs(imag(y))<20eps(T)
    else
        abs(imag(y))/l<20eps(T)
    end
end

issubset(a::Domain,b::Domain) = a==b


first(d::PeriodicDomain) = fromcanonical(d,0)
last(d::PeriodicDomain) = fromcanonical(d,2π)


struct AnyPeriodicDomain <: PeriodicDomain{UnsetNumber} end
isambiguous(::AnyPeriodicDomain)=true

convert(::Type{D},::AnyDomain) where {D<:PeriodicDomain} = AnyPeriodicDomain()

## conveninece routines

ones(d::Domain) = ones(prectype(d),Space(d))
zeros(d::Domain) = zeros(prectype(d),Space(d))





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

commondomain(P::AbstractVector,g::AbstractArray{T}) where {T<:Number} = commondomain(P)
commondomain(P::AbstractVector,g) = commondomain([P;g])


domain(::Number) = AnyDomain()




## rand


rand(d::IntervalDomain,k...) = fromcanonical.(Ref(d),2rand(k...)-1)
rand(d::PeriodicDomain,k...) = fromcanonical.(Ref(d),2π*rand(k...)-π)

checkpoints(d::IntervalDomain) = fromcanonical.(Ref(d),[-0.823972,0.01,0.3273484])
checkpoints(d::PeriodicDomain) = fromcanonical.(Ref(d),[1.223972,3.14,5.83273484])

## boundary

"""
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
invfromcanonicalD(d::Domain,x...) = 1/fromcanonicalD(d,x...)



## domains in higher dimensions


## sorting
# we sort spaces lexigraphically by default

for OP in (:<,:(<=),:>,:(>=),:(isless))
    @eval $OP(a::Domain,b::Domain)=$OP(string(a),string(b))
end


## Other special domains

struct PositiveIntegers <: Domain{Int} end
struct Integers <: Domain{Int} end

const ℕ = PositiveIntegers()
const ℤ = Integers()
