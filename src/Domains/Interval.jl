

export Interval



## Standard interval

immutable Interval{T<:Number} <: IntervalDomain
	a::T
	b::T
end

Interval()=Interval(-1.,1.)
Interval(a::Int,b::Int) = Interval(float64(a),float64(b))   #convenience method

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



Base.convert{D<:Domain}(::Type{D},i::Vector)=Interval(i)
Interval(a::Number,b::Number) = Interval(promote(a,b)...)


## Information

Base.first(d::Interval)=d.a
Base.last(d::Interval)=d.b



## Map interval


#TODO dont know why a union type isn't sufficient to collapse these
tocanonical{T}(d::Interval{T},x::T)=(d.a + d.b - 2x)/(d.a - d.b)
tocanonicalD{T}(d::Interval{T},x::T)=2/( d.b- d.a)
fromcanonical{T}(d::Interval{T},x::T)=(d.a + d.b)/2 + (d.b - d.a)x/2
fromcanonicalD{T}(d::Interval{T},x::T)=( d.b- d.a) / 2

tocanonical{T}(d::Interval{T},x::Array{T})=(d.a + d.b - 2x)/(d.a - d.b)
tocanonicalD{T}(d::Interval{T},x::Array{T})=2/( d.b- d.a)
fromcanonical{T}(d::Interval{T},x::Array{T})=(d.a + d.b)/2 + (d.b - d.a)x/2
fromcanonicalD{T}(d::Interval{T},x::Array{T})=( d.b- d.a) / 2


#possibly a bug prevents me from properly dispatching on x (https://groups.google.com/forum/#!topic/julia-users/_paVdB2wy5k)
for op in (:tocanonical, :tocanonicalD, :fromcanonical, :fromcanonicalD)
    @eval begin 
        function $op{T}(d::Interval{T},x)  
            if typeof(x) <: Array
                return $op(d,convert(Array{T,},x))
            elseif typeof(x) <: Number
                return $op(d,convert(T,x))
            else
                error("Don't know what to do with type: $(typeof(x))")
            end
        end
    end
end


Base.length(d::Interval) = abs(d.b - d.a)
Base.angle(d::Interval)=angle(d.b-d.a)


==(d::Interval,m::Interval) = d.a == m.a && d.b == m.b


##Coefficient space operators

identity_fun(d::Interval)=Fun([.5*(d.b+d.a),.5*(d.b-d.a)],d)


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