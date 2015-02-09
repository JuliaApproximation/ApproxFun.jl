

export Domain,tocanonical,fromcanonical,Circle,PeriodicInterval,Interval
export chebyshevpoints


abstract Domain{T<:Number}  #type parameter represents what find of numeric representation should be used in... TODO explain

immutable AnyDomain <: Domain{UnsetNumber}
end


Base.eltype{T}(::Domain{T})=T



##General routines








## Interval Domains

abstract IntervalDomain{T} <: Domain{T}

function chebyshevpoints{T<:Number}(::Type{T},n::Integer;kind::Integer=1)
    if kind == 1
        return cospi((one(T)/2+[-n:-1])/n)
    elseif kind == 2
        if n==1
            return zeros(T,1)
        else
            return cospi([n-1:-1:0]/(n-one(T)))
        end
    end
end
chebyshevpoints(n::Integer;kind::Integer=1) = chebyshevpoints(Float64,n;kind=kind)

##TODO: Should fromcanonical be fromcanonical!?

points{T}(d::IntervalDomain{T},n::Integer) = fromcanonical(d,chebyshevpoints(T,n))

points(d::Vector,n::Integer)=points(Interval(d),n)
bary(v::Vector{Float64},d::IntervalDomain,x::Float64)=bary(v,tocanonical(d,x))

#TODO consider moving these
Base.first{T}(d::IntervalDomain{T})=fromcanonical(d,-one(T))
Base.last{T}(d::IntervalDomain{T})=fromcanonical(d,one(T))


function Base.in(x,d::IntervalDomain)
    y=tocanonical(d,x)
    abs(imag(y))<10eps() && -1.0-10eps()/length(d)<real(y)<1.0+10eps()/length(d)
end

###### Periodic domains

abstract PeriodicDomain{T} <: Domain{T}

points{T}(d::PeriodicDomain{T},n::Integer) = fromcanonical(d, fourierpoints(n,T))

fourierpoints(n::Integer) = fourierpoints(n,Float64)
fourierpoints{T<:Number}(n::Integer,::Type{T})= convert(T,π)*[-n:2:n-2]/n


function Base.in(x,d::PeriodicDomain)
    y=tocanonical(d,x)
    l=length(d)
    abs(imag(y))/l<20eps() && -π-2*l*eps()<=real(y)<=π+2*l*eps()
end


Base.first(d::PeriodicDomain)=fromcanonical(d,-π)
Base.last(d::PeriodicDomain)=fromcanonical(d,π)

## conveninece routines

Base.ones(d::Domain)=ones(Space(d))
Base.zeros(d::Domain)=zeros(Space(d))





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
