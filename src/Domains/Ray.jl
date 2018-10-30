

export Ray



## Standard interval

# angle is (false==0) and π (true==1)
# or ranges from (-1,1].  We use 1 as 1==true.
#orientation true means oriented out
"""
    Ray{a}(c,o)

represents a ray at angle `a` starting at `c`, with orientation out to
infinity (`o = true`) or back from infinity (`o = false`).
"""
struct Ray{angle,T<:Number} <: SegmentDomain{T}
    center::T
    orientation::Bool
    Ray{angle,T}(c,o) where {angle,T} = new{angle,T}(c,o)
    Ray{angle,T}(c) where {angle,T} = new{angle,T}(c,true)
    Ray{angle,T}() where {angle,T} = new{angle,T}(zero(T),true)
    Ray{angle,T}(r::Ray{angle,T}) where {angle,T} = r
end

const RealRay{T} = Union{Ray{false,T},Ray{true,T}}

Ray{a}(c,o) where {a} = Ray{a,typeof(c)}(c,o)
Ray{a}(c::Number) where {a} = Ray{a,typeof(c)}(c)
Ray{a}() where {a} = Ray{a,Float64}()

angle(d::Ray{a}) where {a} = a*π

# ensure the angle is always in (-1,1]
Ray(c,a,o) = Ray{a==0 ? false : (abs(a)≈(1.0π) ? true : mod(a/π-1,-2)+1),typeof(c)}(c,o)
Ray(c,a) = Ray(c,a,true)

Ray() = Ray{false}()



##deal with vector

function convert(::Type{Ray}, d::AbstractInterval)
    a,b = endpoints(d)
    @assert abs(a)==Inf || abs(b)==Inf

    if abs(b)==Inf
        Ray(a,angle(b),true)
    else #abs(a)==Inf
        Ray(b,angle(a),false)
    end
end
Ray(d::AbstractInterval) = convert(Ray, d)


isambiguous(d::Ray)=isnan(d.center)
convert(::Type{Ray{a,T}},::AnyDomain) where {a,T<:Number} = Ray{a,T}(NaN,true)
convert(::Type{IT},::AnyDomain) where {IT<:Ray} = Ray(NaN,NaN)


isempty(::Ray) = false

## Map interval

function mobiuspars(d::Ray)
    s=(d.orientation ? 1 : -1)
    α=conj(cisangle(d))
    c=d.center
    s*α,-s*(1+α*c),α,1-α*c
end


for OP in (:mobius,:mobiusinv,:mobiusD,:mobiusinvD)
    @eval $OP(a::Ray,z) = $OP(mobiuspars(a)...,z)
end

ray_tocanonical(x) = isinf(x) ? one(x) : (x-1)/(1+x)
ray_tocanonicalD(x) = isinf(x) ? zero(x) : 2*(1/(1+x))^2
ray_fromcanonical(x) = (1+x)/(1-x)
ray_fromcanonicalD(x) = 2*(1/(x-1))^2
ray_invfromcanonicalD(x) = (x-1)^2/2


for op in (:ray_tocanonical,:ray_tocanonicalD)
    @eval $op(o,x)=(o ? 1 : -1)*$op(x)
end
ray_fromcanonical(o,x)=ray_fromcanonical((o ? 1 : -1)*x)
ray_fromcanonicalD(o,x)=(o ? 1 : -1)*ray_fromcanonicalD((o ? 1 : -1)*x)
ray_invfromcanonicalD(o,x)=(o ? 1 : -1)*ray_invfromcanonicalD((o ? 1 : -1)*x)

cisangle(::Ray{a}) where {a}=cis(a*π)
cisangle(::Ray{false})=1
cisangle(::Ray{true})=-1

tocanonical(d::Ray,x) =
    ray_tocanonical(d.orientation,conj(cisangle(d)).*(x-d.center))
tocanonicalD(d::Ray,x) =
    conj(cisangle(d)).*ray_tocanonicalD(d.orientation,conj(cisangle(d)).*(x-d.center))
fromcanonical(d::Ray,x) = cisangle(d)*ray_fromcanonical(d.orientation,x)+d.center
fromcanonical(d::Ray{false},x) = ray_fromcanonical(d.orientation,x)+d.center
fromcanonical(d::Ray{true},x) = -ray_fromcanonical(d.orientation,x)+d.center
fromcanonicalD(d::Ray,x) = cisangle(d)*ray_fromcanonicalD(d.orientation,x)
invfromcanonicalD(d::Ray,x) = conj(cisangle(d))*ray_invfromcanonicalD(d.orientation,x)




arclength(d::Ray) = Inf

==(d::Ray{a},m::Ray{a}) where {a} = d.center == m.center


mappoint(a::Ray{false}, b::Ray{false}, x::Number) =
    x - a.center + b.center


function mappoint(a::Ray, b::Ray, x::Number)
    d = x - a.center;
    k = d * exp((angle(b)-angle(d))*im)
    k + b.center
end
