

export Circle


##  Circle


immutable Circle{T<:Number,V<:Real} <: PeriodicDomain{Complex{V}}
	center::T
	radius::V
	orientation::Bool
end

Circle{T1<:Number,T2<:Number,V<:Real}(::Type{T1},c::T2,r::V,o::Bool) = Circle(convert(promote_type(T1,T2,V),c),
																	  convert(promote_type(real(T1),real(T2),V),r),o)
Circle{T1<:Number}(::Type{T1},c::Number,r::Real)=Circle(T1,c,r)
Circle(c::Number,r::Real)=Circle(c,r,true)
Circle(r::Real) = Circle(zero(r),r)
Circle(r::Int)=Circle(Float64,0.,r)

Circle{V<:Real}(::Type{V}) = Circle(one(V))
Circle()=Circle(1.)



isambiguous(d::Circle)=isnan(d.center) && isnan(d.radius)
Base.convert{T<:Number,V<:Number}(::Type{Circle{T,V}},::AnyDomain)=Circle{T,V}(NaN,NaN)
Base.convert{IT<:Circle}(::Type{IT},::AnyDomain)=Circle(NaN,NaN)


function tocanonical(d::Circle,ζ)
    v=mappoint(d,Circle(),ζ)- 0.im#Subtract 0.im so branch cut is right
    -1.im.*log(v)
end

tocanonicalD(d::Circle,ζ)=-1.im./(ζ-d.center)  #TODO: Check formula
fromcanonical(d::Circle,θ)=d.radius*exp((d.orientation?1:-1)*1.im*θ) + d.center
fromcanonicalD(d::Circle,θ)=(d.orientation?1:-1)*d.radius*1.im*exp((d.orientation?1:-1)*1.im*θ)


Base.in(z,d::Circle)=isapprox(abs(z-d.center),d.radius)

Base.length(d::Circle) = 2π*d.radius
complexlength(d::Circle)=(d.orientation?1:-1)*im*length(d)  #TODO: why?


==(d::Circle,m::Circle) = d.center == m.center && d.radius == m.radius && d.orientation == m.orientation



function mappoint(d1::Circle,d2::Circle,z)
   v=(z-d1.center)/d1.radius
   d1.orientation != d2.orientation && (v=1./v)
   v*d2.radius+d2.center
end



Base.reverse(d::Circle)=Circle(d.center,d.radius,!d.orientation)
Base.conj(d::Circle)=Circle(conj(d.center),d.radius,!d.orientation)


for op in (:+,:-,:.+,:.-)
    @eval begin
        $op(c::Number,d::Circle)=Circle($op(c,d.center),d.radius,d.orientation)
        $op(d::Circle,c::Number)=Circle($op(d.center,c),d.radius,d.orientation)
    end
end

for op in (:*,:.*)
    @eval begin
        $op(c::Real,d::Circle)=Circle($op(c,d.center),$op(abs(c),d.radius),sign(c)<0?!d.orientation:d.orientation)
        $op(d::Circle,c::Real)=Circle($op(c,d.center),$op(abs(c),d.radius),sign(c)<0?!d.orientation:d.orientation)
    end
end
