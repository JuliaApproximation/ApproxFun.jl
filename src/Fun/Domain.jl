

export Domain,tocanonical,fromcanonical,Circle,PeriodicInterval,Interval
export chebyshevpoints


abstract Domain{T<:Number}  #type parameter represents what find of numeric representation should be used in... TODO explain

immutable AnyDomain <: Domain{UnsetNumber}
end


Base.eltype{T}(::Domain{T})=T



##General routines








## Interval Domains

abstract IntervalDomain{T} <: Domain{T}

##TODO: Should fromcanonical be fromcanonical!?
function chebyshevroots(n::Integer,numbertype::Type)
    _π = convert(numbertype,π)
    return [cos(_π*k) for k=-1.+1/(2n):1/n:-1./(2n)]
end

function chebyshevpoints(n::Integer,numbertype::Type)
    if n==1
        return zeros(numbertype,1)
    else
        _π = convert(numbertype,π)
        return [cos(_π*k/(n-1)) for k = n-1:-1:0]
    end
end 

chebyshevpoints(n::Integer) = chebyshevpoints(n,Float64)
chebyshevroots(n::Integer) = chebyshevroots(n,Float64)

function points{T}(d::IntervalDomain{T},n::Integer) 
    if n==1
        return [fromcanonical(d,zero(T))]
    else
        _π = convert(T,π)
        return [fromcanonical(d,cos(_π*k/(n-1))) for k = n-1:-1:0]  #TODO, refactor to use chebyshevpoints
    end
end

points(d::Vector,n::Integer)=points(Interval(d),n)
bary(v::Vector{Float64},d::IntervalDomain,x::Float64)=bary(v,tocanonical(d,x))

#TODO consider moving these
Base.first(d::IntervalDomain)=fromcanonical(d,-1.0)
Base.last(d::IntervalDomain)=fromcanonical(d,1.0)


function Base.in(x,d::IntervalDomain)
    y=tocanonical(d,x)
    abs(imag(y))<10eps() && -1.0-10eps()/length(d)<real(y)<1.0+10eps()/length(d)
end

###### Periodic domains

abstract PeriodicDomain{T} <: Domain{T}

points{T}(d::PeriodicDomain{T},n::Integer) = fromcanonical(d, fourierpoints(n,T))

fourierpoints(n::Integer) = fourierpoints(n,Float64)
fourierpoints(n::Integer,numbertype::Type)= convert(numbertype,π)*[-1.:2/n:1. - 2/n]


function Base.in(x,d::PeriodicDomain)
    y=tocanonical(d,x)
    l=length(d)
    abs(imag(y))/l<20eps() && -π-2*l*eps()<=real(y)<=π+2*l*eps()
end


Base.first(d::PeriodicDomain)=fromcanonical(d,-π)
Base.last(d::PeriodicDomain)=fromcanonical(d,π)

## conveninece routines

Base.ones(d::Domain)=Fun(one,d)
Base.zeros(d::Domain)=Fun(zero,d)





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