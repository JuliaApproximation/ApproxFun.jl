

export Ray



## Standard interval

# angle is (false==0) and π (true==1)
# or ranges from (-1,1].  We use 1 as 1==true.
#orientation true means oriented out
doc"""
    Ray{a}(c,o)

represents a ray at angle `a` starting at `c`, with orientation out to
infinity (`o = true`) or back from infinity (`o = false`).
"""
immutable Ray{angle,T<:Number} <: IntervalDomain{T}
    center::T
    orientation::Bool
    (::Type{Ray{angle,T}}){angle,T}(c,o) = new{angle,T}(c,o)
    (::Type{Ray{angle,T}}){angle,T}(c) = new{angle,T}(c,true)
    (::Type{Ray{angle,T}}){angle,T}() = new{angle,T}(zero(T),true)
    (::Type{Ray{angle,T}}){angle,T}(r::Ray{angle,T}) = r
end

const RealRay{T} = Union{Ray{false,T},Ray{true,T}}

(::Type{Ray{a}}){a}(c,o) = Ray{a,typeof(c)}(c,o)
(::Type{Ray{a}}){a}(c::Number) = Ray{a,typeof(c)}(c)
(::Type{Ray{a}}){a}() = Ray{a,Float64}()

Base.angle{a}(d::Ray{a}) = a*π

# ensure the angle is always in (-1,1]
Ray(c,a,o) = Ray{a==0?false:(abs(a)==(1.0π)?true:mod(a/π-1,-2)+1),typeof(c)}(c,o)
Ray(c,a) = Ray(c,a,true)

Ray() = Ray{false}()



##deal with vector

function Base.convert(::Type{Ray},d::ClosedInterval)
    a,b=d.left,d.right
    @assert abs(a)==Inf || abs(b)==Inf

    if abs(b)==Inf
        Ray(a,angle(b),true)
    else #abs(a)==Inf
        Ray(b,angle(a),false)
    end
end


isambiguous(d::Ray)=isnan(d.center)
Base.convert{a,T<:Number}(::Type{Ray{a,T}},::AnyDomain) = Ray{a,T}(NaN,true)
Base.convert{IT<:Ray}(::Type{IT},::AnyDomain) = Ray(NaN,NaN)




## Map interval

function mobiuspars(d::Ray)
    s=(d.orientation?1:-1)
    α=conj(cisangle(d))
    c=d.center
    s*α,-s*(1+α*c),α,1-α*c
end


for OP in (:mobius,:mobiusinv,:mobiusD,:mobiusinvD)
    @eval $OP(a::Ray,z) = $OP(mobiuspars(a)...,z)
end

ray_tocanonical(x)=(x==Inf)?1.:(x-1.)./(1+x)
ray_tocanonicalD(x)=(x==Inf)?0.:2*(1./(1+x))^2
ray_fromcanonical(x)=(1+x)./(1-x)
ray_fromcanonicalD(x)=2*(1./(x-1.))^2
ray_invfromcanonicalD(x)=(x-1.)^2/2


for op in (:ray_tocanonical,:ray_tocanonicalD)
    @eval $op(o,x)=(o?1:-1)*$op(x)
end
ray_fromcanonical(o,x)=ray_fromcanonical((o?1:-1)*x)
ray_fromcanonicalD(o,x)=(o?1:-1)*ray_fromcanonicalD((o?1:-1)*x)
ray_invfromcanonicalD(o,x)=(o?1:-1)*ray_invfromcanonicalD((o?1:-1)*x)

cisangle{a}(::Ray{a})=cis(a*π)
cisangle(::Ray{false})=1
cisangle(::Ray{true})=-1

tocanonical(d::Ray,x) =
    ray_tocanonical(d.orientation,conj(cisangle(d)).*(x-d.center))
tocanonicalD(d::Ray,x) =
    conj(cisangle(d)).*ray_tocanonicalD(d.orientation,conj(cisangle(d)).*(x-d.center))
fromcanonical(d::Ray,x) = cisangle(d)*ray_fromcanonical(d.orientation,x)+d.center
fromcanonicalD(d::Ray,x) = cisangle(d)*ray_fromcanonicalD(d.orientation,x)
invfromcanonicalD(d::Ray,x) = conj(cisangle(d))*ray_invfromcanonicalD(d.orientation,x)




arclength(d::Ray) = Inf

=={a}(d::Ray{a},m::Ray{a}) = d.center == m.center
