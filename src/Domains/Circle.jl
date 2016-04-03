

export Circle


##  Circle


# T Must be in an Algebra
immutable Circle{T,V<:Real,TT} <: PeriodicDomain{TT}
	center::T
	radius::V
	orientation::Bool

end

Circle(c::Number,r::Real,o::Bool)=Circle{typeof(c),typeof(r),Complex{typeof(r)}}(c,r,o)
Circle(c::Vec,r::Real,o::Bool)=Circle{typeof(c),typeof(r),typeof(c)}(c,r,o)

Circle{T1,T2,V<:Real}(::Type{T1},c::T2,r::V,o::Bool) = Circle(convert(promote_type(T1,T2,V),c),
															  convert(promote_type(real(T1),real(T2),V),r),o)
Circle{T1<:Number}(::Type{T1},c,r::Bool)=Circle(T1,c,r)
Circle{T1<:Number}(::Type{T1},c,r::Real)=Circle(T1,c,r)
Circle(c,r::Real)=Circle(c,r,true)
Circle(r::Real) = Circle(zero(r),r)
Circle(r::Int)=Circle(Float64,0.,r)
Circle(a::Tuple,r::Real)=Circle(Vec(a...),r)

Circle{V<:Real}(::Type{V}) = Circle(one(V))
Circle()=Circle(1.)



isambiguous(d::Circle)=isnan(d.center) && isnan(d.radius)
Base.convert{T<:Number,V<:Real}(::Type{Circle{T,V}},::AnyDomain)=Circle{T,V}(NaN,NaN)
Base.convert{IT<:Circle}(::Type{IT},::AnyDomain)=Circle(NaN,NaN)


function tocanonical{T<:Number}(d::Circle{T},ζ)
    v=mappoint(d,Circle(),ζ)
    atan2(imag(v)-0.0,real(v))  # -0.0 to get branch cut right
end

function tocanonical{T<:Vec}(d::Circle{T},ζ)
    v=mappoint(d,Circle((0.,0.),1.),ζ)
    atan2(v[2]-0.0,v[1])  # -0.0 to get branch cut right
end

fromcanonical{T<:Number}(d::Circle{T},θ)=d.radius*exp((d.orientation?1:-1)*1.im*θ) + d.center
fromcanonicalD{T<:Number}(d::Circle{T},θ)=(d.orientation?1:-1)*d.radius*1.im*exp((d.orientation?1:-1)*1.im*θ)


fromcanonical{T<:Vec}(d::Circle{T},θ::Number)=d.radius*Vec(cos((d.orientation?1:-1)*θ),sin((d.orientation?1:-1)*θ)) + d.center
fromcanonicalD{T<:Vec}(d::Circle{T},θ::Number)=d.radius*(d.orientation?1:-1)*Vec(-sin((d.orientation?1:-1)*θ),cos((d.orientation?1:-1)*θ))


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


/(c::Number,d::Circle)=c==1?(d.center==0?Circle(d.center,1/d.radius,!d.orientation):Circle(1/d.center,abs(1/(d.center+d.radius)-1/(d.center)),!d.orientation)):c*(1/d)

