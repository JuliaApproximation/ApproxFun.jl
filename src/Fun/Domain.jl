

export Domain,tocanonical,fromcanonical,Circle,PeriodicInterval,Interval
export chebyshevpoints


abstract Domain{T<:Number}  #type parameter represents what find of numeric representation should be used in... TODO explain

immutable AnyDomain <: Domain{UnsetNumber} end




Base.eltype{T}(::Domain{T})=T
Base.isreal{T<:Real}(::Domain{T})=true
Base.isreal{T}(::Domain{T})=false

complexlength(::AnyDomain)=NaN
Base.length(::AnyDomain)=NaN


##General routines








## Interval Domains

abstract IntervalDomain{T} <: Domain{T}

canonicaldomain(::IntervalDomain)=Interval()

function chebyshevpoints{T<:Number}(::Type{T},n::Integer;kind::Integer=1)
    if kind == 1
        T[cospi((one(T)/2+k)/n) for k=-n:-1]
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

points{T}(d::IntervalDomain{T},n::Integer) = fromcanonical(d,chebyshevpoints(real(T),n))

points(d::Vector,n::Integer)=points(Interval(d),n)
bary(v::Vector{Float64},d::IntervalDomain,x::Float64)=bary(v,tocanonical(d,x))

#TODO consider moving these
Base.first{T}(d::IntervalDomain{T})=fromcanonical(d,-one(T))
Base.last{T}(d::IntervalDomain{T})=fromcanonical(d,one(T))


function Base.in{T}(x,d::IntervalDomain{T})
    y=tocanonical(d,x)
    abs(imag(y))<100eps(T)/length(d) && -one(real(T))-100eps(T)/length(d)<real(y)<one(real(T))+100eps(T)/length(d)
end

###### Periodic domains

abstract PeriodicDomain{T} <: Domain{T}

points{T}(d::PeriodicDomain{T},n::Integer) = fromcanonical(d, fourierpoints(T,n))

fourierpoints(n::Integer) = fourierpoints(Float64,n)
fourierpoints{T<:Number}(::Type{T},n::Integer)= convert(T,π)*collect(-n:2:n-2)/n


function Base.in{T}(x,d::PeriodicDomain{T})
    y=tocanonical(d,x)
    l=length(d)
    abs(imag(y))/l<20eps(T) && -π-2l*eps(T)<real(y)<π+2l*eps(T)
end


Base.first(d::PeriodicDomain)=fromcanonical(d,-π)
Base.last(d::PeriodicDomain)=fromcanonical(d,π)

## conveninece routines

Base.ones(d::Domain)=ones(eltype(d),Space(d))
Base.zeros(d::Domain)=zeros(eltype(d),Space(d))





function commondomain(P::Vector)
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

commondomain{T<:Number}(P::Vector,g::Array{T})=commondomain(P)
commondomain(P::Vector,g)=commondomain([P,g])


domain(::Number)=AnyDomain()




## rand

Base.rand(d::IntervalDomain,k...)=fromcanonical(d,2rand(k...)-1)
Base.rand(d::PeriodicDomain,k...)=fromcanonical(d,2π*rand(k...)-π)

checkpoints(d::IntervalDomain)=fromcanonical(d,[-0.823972,0.3273484])
checkpoints(d::PeriodicDomain)=fromcanonical(d,[1.223972,-2.83273484])

## boundary


∂(d::IntervalDomain)=[first(d),last(d)]
∂(d::PeriodicDomain)=[]




## map domains


mappoint(d1::Domain,d2::Domain,x)=fromcanonical(d2,tocanonical(d1,x))
