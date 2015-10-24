

export Interval



## Standard interval

immutable Interval{T} <: IntervalDomain{T}  #repeat <:Number due to Julia issue #9441
	a::T
	b::T
	Interval()=new(-one(T),one(T))
	Interval(a,b)=new(a,b)
end

Interval()=Interval{Float64}()
Interval(a::Int,b::Int) = Interval(Float64(a),Float64(b)) #convenience method
Interval(a::Number,b::Number) = Interval{promote_type(typeof(a),typeof(b))}(a,b)

function Interval{T<:Number}(d::AbstractVector{T})
    @assert length(d)==2
    @assert isfinite(d[1]) && isfinite(d[2])
    Interval(d...)
end



Base.convert{T<:Number}(::Type{Interval{T}}, d::Interval) = Interval{T}(d.a,d.b)

AnyInterval{T}(::Type{T})=Interval{T}(NaN,NaN)
AnyInterval()=AnyInterval(Float64)
isambiguous(d::Interval)=isnan(d.a) && isnan(d.b)
Base.convert{T<:Number}(::Type{Interval{T}},::AnyDomain)=AnyInterval(T)
Base.convert{IT<:Interval}(::Type{IT},::AnyDomain)=AnyInterval()


## Information

Base.first(d::Interval)=d.a
Base.last(d::Interval)=d.b
Base.isempty(d::Interval)=isapprox(d.a,d.b;atol=200eps(eltype(d)))

Base.issubset(a::Interval,b::Interval)=first(a)∈b && last(a)∈b



## Map interval


tocanonical(d::Interval,x)=(d.a + d.b - 2x)/(d.a - d.b)
tocanonicalD(d::Interval,x)=2/( d.b- d.a)
fromcanonical(d::Interval,x)=(d.a + d.b)/2 + (d.b - d.a)x/2
fromcanonicalD(d::Interval,x)=( d.b- d.a) / 2


Base.length(d::Interval) = abs(d.b - d.a)
Base.angle(d::Interval)=angle(d.b-d.a)
complexlength(d::Interval)=d.b-d.a


==(d::Interval,m::Interval) = d.a == m.a && d.b == m.b
function Base.isapprox(d::Interval,m::Interval)
    tol=10E-12
    abs(d.a-m.a)<tol&&abs(d.b-m.b)<tol
end

##Coefficient space operators

identity_fun(d::Interval)=Fun(eltype(d)[(d.b+d.a)/2,(d.b-d.a)/2],d)



## algebra

for op in (:*,:+,:-,:.*,:.+,:.-,:.^)
    @eval begin
        $op(c::Number,d::Interval)=Interval($op(c,d.a),$op(c,d.b))
        $op(d::Interval,c::Number)=Interval($op(d.a,c),$op(d.b,c))
    end
end

for op in (:/,:./)
    @eval $op(d::Interval,c::Number)=Interval($op(d.a,c),$op(d.b,c))
end

Base.sqrt(d::Interval)=Interval(sqrt(d.a),sqrt(d.b))

+(d1::Interval,d2::Interval)=Interval(d1.a+d2.a,d1.b+d2.b)



## intersect/union

Base.reverse(d::Interval)=Interval(d.b,d.a)

function Base.intersect{T<:Real,V<:Real}(a::Interval{T},b::Interval{V})
    if first(a) > last(a)
        intersect(reverse(a),b)
    elseif first(b) > last(b)
        intersect(a,reverse(b))
    elseif first(a) > first(b)
        intersect(b,a)
    elseif last(a) <= first(b)
        EmptyDomain()
    elseif last(a)>=last(b)
        b
    elseif isapprox(first(b),last(a);atol=100eps(promote_type(T,V))/max(length(a),length(b)))
        EmptyDomain()
    else
        Interval(first(b),last(a))
    end
end


function Base.setdiff{T<:Real,V<:Real}(a::Interval{T},b::Interval{V})
    # ensure a/b are well-ordered
    if first(a) > last(a)
        intersect(reverse(a),b)
    elseif first(b) > last(b)
        intersect(a,reverse(b))
    elseif first(a)< first(b)
        if last(a) <= first(b)
            a
        else # first(a) ≤ first(b) ≤last(a)
            #TODO: setdiff in the middle
            @assert last(a) <= last(b)
            Interval(first(a),first(b))
        end
    else #first(a)>= first(b)
        if first(a)>=last(b)
            a
        elseif last(a) <= last(b)
            EmptyDomain()
        else #first(b) < first(a) < last(b) < last(a)
            Interval(last(b),last(a))
        end
    end
end






## Multivariate

tocanonical{V<:Vector}(d::Interval{V},x)=first(d.a + d.b - 2x)/first(d.a - d.b)
fromcanonical{V<:Vector}(d::Interval{V},p::AbstractVector)=V[fromcanonical(d,x) for x in p]
