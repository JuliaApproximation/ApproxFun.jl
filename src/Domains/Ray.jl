

export Ray



## Standard interval

# TODO: allow bool for angle a la line
immutable Ray{CT<:Number,T<:Number} <: IntervalDomain{T}
    center::CT
    angle::Float64
    orientation::Bool
    Ray(c,a,o)=new(c,a,o)
end


Ray{T<:Number}(c::T,a::Real,o::Bool)=(a==0||a==π)?Ray{T,T}(c,a,o):Ray{T,promote_type(T,Complex{Float64})}(c,a,o)
Ray(c,a,o::Int)=Ray(c,a,o==1)
Ray(c,a)=Ray(c,a,true)
Ray()=Ray(0.,0.)

##deal with vector

function Ray(d::Vector)
    @assert length(d)==2
    @assert abs(d[1])==Inf|| abs(d[2])==Inf

    if abs(d[2])==Inf
        Ray(d[1],angle(d[2]),true)
    else #d[1]==Inf
        Ray(d[2],angle(d[1]),false)
    end
end


isambiguous(d::Ray)=isnan(d.center) && isnan(d.angle)
Base.convert{CT<:Number,T<:Number}(::Type{Ray{CT,T}},::AnyDomain)=Ray{CT,T}(NaN,NaN,true)
Base.convert{IT<:Ray}(::Type{IT},::AnyDomain)=Ray(NaN,NaN)




## Map interval



ray_tocanonical(x)=(x==Inf)?1.:(x-1.)./(1+x)
ray_tocanonicalD(x)=(x==Inf)?0.:2*(1./(1+x)).^2
ray_fromcanonical(x)=(1+x)./(1-x)
ray_fromcanonicalD(x)=2*(1./(x-1.)).^2
ray_invfromcanonicalD(x)=(x-1.).^2/2

# atomatically vectorize over vector arg
@vectorize_1arg Number ray_tocanonical
@vectorize_1arg Number ray_tocanonicalD
@vectorize_1arg Number ray_fromcanonical
@vectorize_1arg Number ray_fromcanonicalD


for op in (:ray_tocanonical,:ray_tocanonicalD)
    @eval $op(o,x)=(o?1:-1)*$op(x)
end
ray_fromcanonical(o,x)=ray_fromcanonical((o?1:-1)*x)
ray_fromcanonicalD(o,x)=(o?1:-1)*ray_fromcanonicalD((o?1:-1)*x)
ray_invfromcanonicalD(o,x)=(o?1:-1)*ray_invfromcanonicalD((o?1:-1)*x)

function cistyped(a::Real)
    if a==0
        1.
    elseif a==π||a==-π
        -1.
    else
        cis(a)
    end
end

tocanonical(d::Ray,x)=ray_tocanonical(d.orientation,cistyped(-d.angle).*(x-d.center))
tocanonicalD(d::Ray,x)=cistyped(-d.angle).*ray_tocanonicalD(d.orientation,cistyped(-d.angle).*(x-d.center))
fromcanonical(d::Ray,x)=cistyped(d.angle)*ray_fromcanonical(d.orientation,x)+d.center
fromcanonicalD(d::Ray,x)=cistyped(d.angle)*ray_fromcanonicalD(d.orientation,x)
invfromcanonicalD(d::Ray,x)=cistyped(-d.angle)*ray_invfromcanonicalD(d.orientation,x)




Base.length(d::Ray) = Inf



==(d::Ray,m::Ray) = d.center == m.center && d.angle == m.angle





