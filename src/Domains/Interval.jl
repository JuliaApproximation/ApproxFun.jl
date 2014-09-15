

export Interval



## Standard interval

immutable Interval{T<:Number} <: IntervalDomain
	a::T
	b::T
end

Interval()=Interval(-1.,1.)


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
        [Interval(d[1:2]),Interval(d[2:end])]
    end
end



Base.convert{D<:Interval}(::Type{D},i::Vector)=Interval(i)
Interval(a::Number,b::Number) = Interval(promote(a,b)...)


## Information

Base.first(d::Interval)=d.a
Base.last(d::Interval)=d.b



## Map interval



tocanonical(d::Interval,x)=(d.a + d.b .- 2x)/(d.a - d.b)
tocanonicalD(d::Interval,x)=2/( d.b- d.a)
fromcanonical(d::Interval,x)=.5*(d.a + d.b) .+ .5*(d.b - d.a)x
fromcanonicalD(d::Interval,x)=.5*( d.b- d.a)



Base.length(d::Interval) = d.b - d.a



==(d::Interval,m::Interval) = d.a == m.a && d.b == m.b


##Coefficient space operators

identity_fun(d::Interval)=Fun([.5*(d.b+d.a),.5*(d.b-d.a)],d)


# function multiplybyx{T<:Number,D<:Interval}(f::IFun{T,UltrasphericalSpace{D}})
#     a = domain(f).a
#     b = domain(f).b
#     g = IFun([0,1,.5*ones(length(f)-1)].*[0,f.coefficients]+[.5*f.coefficients[2:end],0,0],f.space) #Gives multiplybyx on unit interval
#     (b-a)/2*g + (b+a)/2
# end
