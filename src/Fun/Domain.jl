

export Domain,IntervalDomain,PeriodicDomain,tocanonical,fromcanonical,fromcanonicalD,∂
export chebyshevpoints,fourierpoints,isambiguous


# T is the numeric type used to represent the domain
# d is the dimension
abstract Domain{T,d}
typealias UnivariateDomain{T} Domain{T,1}
typealias BivariateDomain{T} Domain{T,2}


Base.eltype{T}(::Domain{T})=T
Base.eltype{T,d}(::Type{Domain{T,d}})=T
Base.isreal{T<:Real}(::Domain{T})=true
Base.isreal{T}(::Domain{T})=false
Base.ndims{T,d}(::Domain{T,d})=d
Base.ndims{T,d}(::Type{Domain{T,d}})=d
Base.ndims{DT<:Domain}(::Type{DT})=ndims(super(DT))


#TODO: bivariate AnyDomain
immutable AnyDomain <: Domain{UnsetNumber} end
immutable EmptyDomain <: Domain{UnsetNumber} end

isambiguous(::AnyDomain)=true
Base.ndims(::AnyDomain)=1

complexlength(::AnyDomain)=NaN
Base.length(::AnyDomain)=NaN

Base.reverse(a::Union{AnyDomain,EmptyDomain})=a

canonicaldomain(a::Union{AnyDomain,EmptyDomain})=a

Base.in(x::Domain,::EmptyDomain)=false

##General routines


Base.isempty(::EmptyDomain)=true
Base.isempty(::Domain)=false
Base.intersect(::Domain,::Domain)=EmptyDomain()
\(a::Domain,b::Domain)=setdiff(a,b)

## Interval Domains

abstract IntervalDomain{T} <: UnivariateDomain{T}

canonicaldomain(::IntervalDomain)=Interval()

Base.isapprox(a::Domain,b::Domain)=a==b
domainscompatible(a,b) = domainscompatible(domain(a),domain(b))
domainscompatible(a::Domain,b::Domain)=isambiguous(a) || isambiguous(b) ||
                    isapprox(a,b) || isapprox(a,reverse(b))

function chebyshevpoints{T<:Number}(::Type{T},n::Integer;kind::Integer=1)
    if kind == 1
        T[sinpi((n-2k-one(T))/2n) for k=n-1:-1:0]
    elseif kind == 2
        if n==1
            zeros(T,1)
        else
            T[cospi(k/(n-one(T))) for k=n-1:-1:0]
        end
    end
end
chebyshevpoints(n::Integer;kind::Integer=1) = chebyshevpoints(Float64,n;kind=kind)

##TODO: Should fromcanonical be fromcanonical!?

points{T}(d::IntervalDomain{T},n::Integer) = fromcanonical(d,chebyshevpoints(real(eltype(T)),n))  # eltype to handle point

points(d::AbstractVector,n::Integer)=points(convert(Domain,d),n)
bary(v::AbstractVector{Float64},d::IntervalDomain,x::Float64)=bary(v,tocanonical(d,x))

#TODO consider moving these
Base.first{T}(d::IntervalDomain{T})=fromcanonical(d,-one(T))
Base.last{T}(d::IntervalDomain{T})=fromcanonical(d,one(T))

Base.in(x,::AnyDomain)=true
function Base.in{T}(x,d::IntervalDomain{T})
    y=tocanonical(d,x)
    ry=real(y)
    sc=abs(fromcanonicalD(d,ry<-1?-1:(ry>1?1:ry)))  # scale based on stretch of map on projection to interal
    abs(imag(y))≤100eps(T)/sc && -one(real(T))-100eps(T)/sc≤ry≤one(real(T))+100eps(T)/sc
end

pieces(d::Domain)=[d]
issubcomponent(a::Domain,b::Domain)=a in pieces(b)

###### Periodic domains

abstract PeriodicDomain{T} <: UnivariateDomain{T}


canonicaldomain(::PeriodicDomain)=PeriodicInterval()

points{T}(d::PeriodicDomain{T},n::Integer) = fromcanonical(d, fourierpoints(T,n))

fourierpoints(n::Integer) = fourierpoints(Float64,n)
fourierpoints{T<:Number}(::Type{T},n::Integer)= convert(T,π)*collect(-n:2:n-2)/n


function Base.in{T}(x,d::PeriodicDomain{T})
    y=tocanonical(d,x)
    l=length(d)
    if isinf(l)
        abs(imag(y))<20eps(T) && -π-2eps(T)<real(y)<π+2eps(T)
    else
        abs(imag(y))/l<20eps(T) && -π-2l*eps(T)<real(y)<π+2l*eps(T)
    end
end

Base.issubset(a::Domain,b::Domain)=a==b


Base.first(d::PeriodicDomain)=fromcanonical(d,-π)
Base.last(d::PeriodicDomain)=fromcanonical(d,π)


immutable AnyPeriodicDomain <: PeriodicDomain{UnsetNumber} end
isambiguous(::AnyPeriodicDomain)=true

Base.convert{D<:PeriodicDomain}(::Type{D},::AnyDomain)=AnyPeriodicDomain()

## conveninece routines

Base.ones(d::Domain)=ones(eltype(d),Space(d))
Base.zeros(d::Domain)=zeros(eltype(d),Space(d))





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

commondomain{T<:Number}(P::AbstractVector,g::Array{T})=commondomain(P)
commondomain(P::AbstractVector,g)=commondomain([P;g])


domain(::Number)=AnyDomain()




## rand

Base.rand(d::IntervalDomain,k...)=fromcanonical(d,2rand(k...)-1)
Base.rand(d::PeriodicDomain,k...)=fromcanonical(d,2π*rand(k...)-π)

checkpoints(d::IntervalDomain)=fromcanonical(d,eltype(d)[-0.823972,0.3273484])
checkpoints(d::PeriodicDomain)=fromcanonical(d,eltype(d)[1.223972,-2.83273484])

## boundary

∂(d::Domain)=EmptyDomain()   # This is meant to be overriden
∂(d::IntervalDomain)=[first(d),last(d)] #TODO: Points domain
∂(d::PeriodicDomain)=EmptyDomain()




## map domains

fromcanonical(d::Domain,v::AbstractMatrix)=[fromcanonical(d,vk) for vk in v]

mappoint(d1::Domain,d2::Domain,x...)=fromcanonical(d2,tocanonical(d1,x...))
invfromcanonicalD(d::Domain,x...)=1./fromcanonicalD(d,x...)



## domains in higher dimensions

points{T<:Array}(d::IntervalDomain{T},n::Integer) = T[fromcanonical(d,x) for x in chebyshevpoints(real(eltype(T)),n)]


## sorting
# we sort spaces lexigraphically by default

for OP in (:<,:(<=),:>,:(>=),:(Base.isless))
    @eval $OP(a::Domain,b::Domain)=$OP(string(a),string(b))
end
