

export Interval



## Standard interval

immutable Interval{T<:Number} <: IntervalDomain{T}  #repeat <:Number due to Julia issue #9441
	a::T
	b::T
	Interval()=new(-one(T),one(T))
	Interval(a,b)=new(a,b)
end

Interval()=Interval{Float64}()
Interval{T}(a::T,b::T)=Interval{T}(a,b)
Interval(a::Int,b::Int) = Interval(@compat(Float64(a)),@compat(Float64(b)))   #convenience method

#TODO: in 0.4 change to Domain(d::Vector{T})
function Interval{T<:Number}(d::Vector{T})
    @assert length(d) >1

    if length(d) == 2
        if abs(d[1]) == Inf && abs(d[2]) == Inf
            Line(d)
        elseif abs(d[2]) == Inf || abs(d[1]) == Inf
            Ray(d)
        else
            Interval(d[1],d[2])
        end
    else
        [Interval(d[1:2]);Interval(d[2:end])]   #TODO ensure all Intervals are of the same type
    end
end

function Interval{T1<:Number,T2<:Number}(a::Vector{T1},b::Vector{T2})
    @assert length(a) == length(b)
    if length(a) == 1
        Interval(a[1],b[1])
    else
        [Interval(a[1],b[1]);Interval(a[2:end],b[2:end])]
    end
end

Base.convert{T<:Number}(::Type{Interval{T}}, d::Interval) = Interval{T}(d.a,d.b)

# These are needed for spaces to auto-convert [a,b] to Interval
Base.convert(::Type{AnyDomain},i::Vector)=Interval(i)
Base.convert{D<:IntervalDomain}(::Type{D},i::Vector)=Interval(i)
Base.convert(::Type{Union(IntervalDomain,AnyDomain)},i::Vector)=Interval(i)
Base.convert(::Type{Union(Interval,AnyDomain)},i::Vector)=Interval(i)
Interval(a::Number,b::Number) = Interval{promote_type(typeof(a),typeof(b))}(a,b)

AnyInterval{T}(::Type{T})=Interval{T}(NaN,NaN)
AnyInterval()=AnyInterval(Float64)
isambiguous(d::Interval)=isnan(d.a) && isnan(d.b)
Base.convert{T<:Number}(::Type{Interval{T}},::AnyDomain)=AnyInterval(T)
Base.convert{IT<:Interval}(::Type{IT},::AnyDomain)=AnyInterval()


## Information

Base.first(d::Interval)=d.a
Base.last(d::Interval)=d.b
Base.isempty(d::Interval)=isapprox(d.a,d.b)


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


# function multiplybyx{T<:Number,D<:Interval}(f::IFun{T,UltrasphericalSpace{D}})
#     a = domain(f).a
#     b = domain(f).b
#     g = IFun([0,1,.5*ones(length(f)-1)].*[0,f.coefficients]+[.5*f.coefficients[2:end],0,0],f.space) #Gives multiplybyx on unit interval
#     (b-a)/2*g + (b+a)/2
# end



## algebra

for op in (:*,:+,:-,:.*,:.+,:.-)
    @eval begin
        $op(c::Number,d::Interval)=Interval($op(c,d.a),$op(c,d.b))
        $op(d::Interval,c::Number)=Interval($op(d.a,c),$op(d.b,c))
    end
end

for op in (:/,:./)
    @eval $op(d::Interval,c::Number)=Interval($op(d.a,c),$op(d.b,c))
end


+(d1::Interval,d2::Interval)=Interval(d1.a+d2.a,d1.b+d2.b)



## intersect/union

Base.reverse(d::Interval)=Interval(d.b,d.a)

function Base.intersect(a::Interval{Float64},b::Interval{Float64})
    if first(a) > last(a)
        intersect(reverse(a),b)
    elseif first(b) > last(b)
        intersect(a,reverse(b))
    elseif first(a) > first(b)
        intersect(b,a)
    elseif last(a) <= first(b)
        []
    elseif last(a)>=last(b)
        b
    else
        Interval(first(b),last(a))
    end
end


function Base.setdiff(a::Interval{Float64},b::Interval{Float64})
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
            []
        else #first(b) < first(a) < last(b) < last(a)
            Interval(last(b),last(a))
        end
    end
end

# function Base.sort(d::Vector{Interval{Float64}})
#
# end
