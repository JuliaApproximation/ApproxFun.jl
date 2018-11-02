

export Circle


##  Circle


# T Must be in an Algebra
"""
    Circle(c,r,o)

represents the circle centred at `c` with radius `r` which is positively (`o=true`)
or negatively (`o=false`) oriented.
"""
struct Circle{T,V<:Real,TT} <: PeriodicDomain{TT}
	center::T
	radius::V
	orientation::Bool

end

Circle(c::Number,r::Real,o::Bool) = Circle{typeof(c),typeof(r),Complex{typeof(r)}}(c,r,o)
Circle(c::Vec,r::Real,o::Bool) = Circle{typeof(c),typeof(r),typeof(c)}(c,r,o)

Circle(::Type{T1},c::T2,r::V,o::Bool) where {T1,T2,V<:Real} = Circle(convert(promote_type(T1,T2,V),c),
															  convert(promote_type(real(T1),real(T2),V),r),o)
Circle(::Type{T1},c,r::Bool) where {T1<:Number} = Circle(T1,c,r)
Circle(::Type{T1},c,r::Real) where {T1<:Number} = Circle(T1,c,r)
Circle(c,r::Real) = Circle(c,r,true)
Circle(r::Real) = Circle(zero(r),r)
Circle(r::Int) = Circle(Float64,0.,r)
Circle(a::Tuple,r::Real) = Circle(Vec(a...),r)

Circle(::Type{V}) where {V<:Real} = Circle(one(V))
Circle() = Circle(1.0)



isambiguous(d::Circle{T}) where {T<:Number} = isnan(d.center) && isnan(d.radius)
isambiguous(d::Circle{T}) where {T<:Vec} = all(isnan,d.center) && isnan(d.radius)
convert(::Type{Circle{T,V}},::AnyDomain) where {T<:Number,V<:Real} = Circle{T,V}(NaN,NaN)
convert(::Type{IT},::AnyDomain) where {IT<:Circle} = Circle(NaN,NaN)


function tocanonical(d::Circle{T},ζ) where T<:Number
    v=mappoint(d,Circle(),ζ)
    atan(imag(v)-0.0,real(v))  # -0.0 to get branch cut right
end

function tocanonical(d::Circle{T},ζ) where T<:Vec
    v=mappoint(d,Circle((0.0,0.0),1.0),ζ)
    atan(v[2]-0.0,v[1])  # -0.0 to get branch cut right
end


fromcanonical(d::Circle{T,V,Complex{V}},θ) where {T<:Number,V<:Real} =
	d.radius*exp((d.orientation ? 1 : -1)*1.0im*θ) + d.center
fromcanonicalD(d::Circle{T},θ) where {T<:Number} =
	(d.orientation ? 1 : -1)*d.radius*1.0im*exp((d.orientation ? 1 : -1)*1.0im*θ)


fromcanonical(d::Circle{T},θ::Number) where {T<:Vec} =
	d.radius*Vec(cos((d.orientation ? 1 : -1)*θ),sin((d.orientation ? 1 : -1)*θ)) + d.center
fromcanonicalD(d::Circle{T},θ::Number) where {T<:Vec} =
	d.radius*(d.orientation ? 1 : -1)*Vec(-sin((d.orientation ? 1 : -1)*θ),cos((d.orientation ? 1 : -1)*θ))


indomain(z,d::Circle) = norm(z-d.center) ≈ d.radius

arclength(d::Circle) = 2π*d.radius
complexlength(d::Circle) = (d.orientation ? 1 : -1)*im*arclength(d)  #TODO: why?


==(d::Circle,m::Circle) = d.center == m.center && d.radius == m.radius && d.orientation == m.orientation



mappoint(d1::Circle{T},d2::Circle{V},z) where {T<:Vec,V<:Number} =
	mappoint(Circle(complex(d1.center...),d1.radius),d2,z[1]+im*z[2])

mappoint(d1::Circle{T},d2::Circle{V},z) where {T<:Number,V<:Vec} =
	mappoint(Circle(Vec(d1.center...),d1.radius),d2,Vec(real(z),imag(z)))

function mappoint(d1::Circle,d2::Circle,z)
   v=(z-d1.center)/d1.radius
   d1.orientation != d2.orientation && (v=1/v)
   v*d2.radius+d2.center
end



reverseorientation(d::Circle) = Circle(d.center,d.radius,!d.orientation)
conj(d::Circle) = Circle(conj(d.center),d.radius,!d.orientation)


for op in (:+,:-)
    @eval begin
        $op(c::Number,d::Circle) = Circle($op(c,d.center),d.radius,d.orientation)
        $op(d::Circle,c::Number) = Circle($op(d.center,c),d.radius,d.orientation)
    end
end


*(c::Real,d::Circle) = Circle(*(c,d.center),*(abs(c),d.radius),sign(c)<0 ? !d.orientation : d.orientation)
*(d::Circle,c::Real) = Circle(*(c,d.center),*(abs(c),d.radius),sign(c)<0 ? !d.orientation : d.orientation)



/(c::Number,d::Circle) =
	c==1 ? (d.center==0 ? Circle(d.center,1/d.radius,!d.orientation) :
				Circle(1/d.center,abs(1/(d.center+d.radius)-1/(d.center)),!d.orientation)) :
				c*(1/d)
