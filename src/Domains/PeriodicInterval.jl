

export PeriodicInterval


type PeriodicInterval{T<:Number} <: PeriodicDomain
	a::T
	b::T
end

PeriodicInterval()=PeriodicInterval(-1.π,1.π)


Interval(d::PeriodicInterval)=Interval(d.a,d.b)
PeriodicInterval(d::Interval)=PeriodicInterval(d.a,d.b)




## Map periodic interval


tocanonical(d::PeriodicInterval,x)=1.π.*tocanonical(Interval(d),x)
tocanonicalD(d::PeriodicInterval,x)=1.π.*tocanonicalD(Interval(d),x)
fromcanonical(d::PeriodicInterval,θ)=fromcanonical(Interval(d),θ/π)
fromcanonicalD(d::PeriodicInterval,θ)=fromcanonicalD(Interval(d),θ/π)/π



Base.length(d::PeriodicInterval) = d.b - d.a



==(d::PeriodicInterval,m::PeriodicInterval) = d.a == m.a && d.b == m.b



##Differentiation and integration


function Base.diff{T<:Number,M<:PeriodicInterval}(f::FFun{T,M}) 
    tocanonicalD(f.domain,0)*FFun(
                    ShiftVector(1.im*[firstindex(f.coefficients):-1],
                                1.im*[0:lastindex(f.coefficients)]).*f.coefficients,
                    f.domain)
end

function Base.diff(f::FFun,k::Integer)
    @assert k >= 0
    (k==0)?f:diff(diff(f),k-1)
end




function integrate{T<:Number,M<:PeriodicInterval}(f::FFun{T,M}) 
    tol = 10eps()
    @assert abs(f.coefficients[0]) < tol
    
    ##TODO: mapped domains
    
    @assert f.domain.a ==-π
    @assert f.domain.b ==π        
    FFun(
                    ShiftVector(-1.im./[firstindex(f.coefficients):-1],
                                [0,(-1.im./[1:lastindex(f.coefficients)])]).*f.coefficients,
                    f.domain)
end

Base.sum{T<:Number,M<:PeriodicInterval}(f::FFun{T,M})=f.coefficients[0].*length(f.domain)



