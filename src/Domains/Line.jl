

export Line, PeriodicLine




## Standard interval



# angle is π*a where a is (false==0) and (true==1)
# or ranges from (-1,1].  We use 1 as 1==true.

"""
    Line{a}(c)

represents the line at angle `a` in the complex plane, centred at `c`.
"""
struct Line{angle,T<:Number} <: SegmentDomain{T}
    center::T
    α::Float64
    β::Float64

    #TODO get this inner constructor working again.
    Line{angle,T}(c,α,β) where {angle,T} = new{angle,T}(c,α,β)
    Line{angle,T}(c) where {angle,T} = new{angle,T}(c,-1.,-1.)
    Line{angle,T}() where {angle,T} = new{angle,T}(zero(T),-1.,-1.)
end

const RealLine{T} = Union{Line{false,T},Line{true,T}}

Line{a}(c,α,β) where {a} = Line{a,typeof(c)}(c,α,β)
Line{a}(c::Number) where {a} = Line{a,typeof(c)}(c)
Line{a}() where {a} = Line{a,Float64}()

angle(d::Line{a}) where {a} = a*π

reverseorientation(d::Line{true}) = Line{false}(d.center,d.β,d.α)
reverseorientation(d::Line{false}) = Line{true}(d.center,d.β,d.α)
reverseorientation(d::Line{a}) where {a} = Line{a-1}(d.center,d.β,d.α)

# ensure the angle is always in (-1,1]
Line(c,a,α,β) = Line{mod(a/π-1,-2)+1,typeof(c)}(c,α,β)
Line(c,a) = Line(c,a,-1.,-1.)
# true is negative orientation, false is positive orientation
# this is because false==0 and we take angle=0
Line(b::Bool) = Line{b}()
Line() = Line(false)


isambiguous(d::Line)=isnan(d.center)
convert(::Type{Domain{T}},d::Line{a}) where {a,T<:Number} = Line{a,T}(d.center, d.α, d.β)
convert(::Type{Line{a,T}},::AnyDomain) where {a,T<:Number} = Line{a,T}(NaN)
convert(::Type{IT},::AnyDomain) where {IT<:Line}=Line(NaN,NaN)

## Map interval


##TODO non-1 alpha,beta


function line_tocanonical(α,β,x)
    @assert α==β==-1. || α==β==-.5

    if α==β==-1.
        2x/(1+sqrt(1+4x^2))
    elseif α==β==-.5
        x/sqrt(1 + x^2)
    end
end

function line_tocanonicalD(α,β,x)
    @assert α==β==-1. || α==β==-.5

    if α==β==-1.
        2/(1+4x^2+sqrt(1+4x^2))
    elseif α==β==-0.5
        (1 + x^2)^(-3/2)
    end
end
function line_fromcanonical(α,β,x)
    #TODO: why is this consistent?
    if α==β==-1.
        x/(1-x^2)
    else
        x*(1 + x)^α*(1 - x)^β
    end
end
function line_fromcanonicalD(α,β,x)
    if α==β==-1.
        (1+x^2)/(1-x^2)^2
    else
        (1 - (β-α)x - (β+α+1)x^2)*(1+x)^(α-1)*(1-x)^(β-1)
    end
end

function line_invfromcanonicalD(α,β,x)
    if α==β==-1.
        (1-x^2)^2/(1+x^2)
    else
        1/(1 - (β-α)x - (β+α+1)x^2)*(1+x)^(1-α)*(1-x)^(1-β)
    end
end


tocanonical(d::Line,x) = line_tocanonical(d.α,d.β,cis(-angle(d)).*(x-d.center))
tocanonical(d::Line{false},x) = line_tocanonical(d.α,d.β,x-d.center)
tocanonical(d::Line{true},x) = line_tocanonical(d.α,d.β,d.center-x)

tocanonicalD(d::Line,x) = cis(-angle(d)).*line_tocanonicalD(d.α,d.β,cis(-angle(d)).*(x-d.center))
tocanonicalD(d::Line{false},x) = line_tocanonicalD(d.α,d.β,x-d.center)
tocanonicalD(d::Line{true},x) = -line_tocanonicalD(d.α,d.β,d.center-x)

fromcanonical(d::Line,x) = cis(angle(d))*line_fromcanonical(d.α,d.β,x)+d.center
fromcanonical(d::Line{false},x) = line_fromcanonical(d.α,d.β,x)+d.center
fromcanonical(d::Line{true},x) = -line_fromcanonical(d.α,d.β,x)+d.center

fromcanonicalD(d::Line,x) = cis(angle(d))*line_fromcanonicalD(d.α,d.β,x)
fromcanonicalD(d::Line{false},x) = line_fromcanonicalD(d.α,d.β,x)
fromcanonicalD(d::Line{true},x) = -line_fromcanonicalD(d.α,d.β,x)

invfromcanonicalD(d::Line,x) = cis(-angle(d))*line_invfromcanonicalD(d.α,d.β,x)
invfromcanonicalD(d::Line{false},x) = line_invfromcanonicalD(d.α,d.β,x)
invfromcanonicalD(d::Line{true},x) = -line_invfromcanonicalD(d.α,d.β,x)






==(d::Line{a},m::Line{a}) where {a} = d.center == m.center && d.β == m.β &&d.α == m.α



# algebra
*(c::Real,d::Line{false}) = Line{sign(c)>0 ? false : true}(isapprox(d.center,0) ? d.center : c*d.center,d.α,d.β)
*(c::Real,d::Line{true}) = Line{sign(c)>0 ? true : false}(isapprox(d.center,0) ? d.center : c*d.center,d.α,d.β)
*(c::Number,d::Line) = Line(isapprox(d.center,0) ? d.center : c*d.center,angle(d)+angle(c),d.α,d.β)
*(d::Line,c::Number) = c*d
for OP in (:+,:-)
    @eval begin
        $OP(c::Number,d::Line{a}) where {a} = Line{a}($OP(c,d.center),d.α,d.β)
        $OP(d::Line{a},c::Number) where {a} = Line{a}($OP(d.center,c),d.α,d.β)
    end
end






## Periodic line

# angle is (false==0) and π (true==1)
# or ranges from (-1,1]
struct PeriodicLine{angle,T} <: PeriodicDomain{T}
    center::T
    L::Float64
    PeriodicLine{angle,T}(c,L) where {angle,T} = new{angle,T}(c,L)
    PeriodicLine{angle,T}(c) where {angle,T} = new{angle,T}(c,1.)
    PeriodicLine{angle,T}(d::PeriodicLine) where {angle,T} = new{angle,T}(d.center,d.L)
    PeriodicLine{angle,T}() where {angle,T} = new{angle,T}(0.,1.)
end

(::Type{PeriodicLine{a}})(c,L) where {a} = PeriodicLine{a,typeof(c)}(c,L)


PeriodicLine(c,a) = PeriodicLine{a/π,eltype(c)}(c,1.)
PeriodicLine() = PeriodicLine{false,Float64}(0.,1.)
PeriodicLine(b::Bool) = PeriodicLine{b,Float64}()

isambiguous(d::PeriodicLine) = isnan(d.center) && isnan(d.angle)
convert(::Type{Domain{T}},d::PeriodicLine{a}) where {a,T<:Number} = PeriodicLine{a,T}(d.center, d.L)
convert(::Type{PeriodicLine{T,TT}},::AnyDomain) where {T<:Number,TT} = PeriodicLine{T,TT}(NaN,NaN)
convert(::Type{IT},::AnyDomain) where {IT<:PeriodicLine} = PeriodicLine(NaN,NaN)

angle(d::PeriodicLine{a}) where {a}=a*π

reverseorientation(d::PeriodicLine{true})=PeriodicLine{false}(d.center,d.L)
reverseorientation(d::PeriodicLine{false})=PeriodicLine{true}(d.center,d.L)
reverseorientation(d::PeriodicLine{a}) where {a}=PeriodicLine{a-1}(d.center,d.L)

tocanonical(d::PeriodicLine{false},x) = real(2atan((x-d.center)/d.L))
fromcanonical(d::PeriodicLine{false},θ) = d.L*tan(θ/2) + d.center

tocanonical(d::PeriodicLine{a},x) where {a} = tocanonical(PeriodicLine{false,Float64}(0.,d.L),exp(-π*im*a)*(x-d.center))
fromcanonical(d::PeriodicLine{a},v::AbstractArray) where {a} =
    [fromcanonical(d,vk) for vk in v]
fromcanonical(d::PeriodicLine{a},x) where {a} =
    exp(π*im*a)*fromcanonical(PeriodicLine{false,Float64}(0.,d.L),x)+d.center


function invfromcanonicalD(d::PeriodicLine{false})
    @assert d.center==0  && d.L==1.0
    a=Fun(PeriodicInterval(),[1.,0,1])
end

mappoint(a::PeriodicLine{false},b::Circle,x) = b.radius*((a.L*im-(x-a.center))./(a.L*im+(x-a.center)))+b.center
function mappoint(b::Circle,a::PeriodicLine{false},x)
    y=(x-b.center)./b.radius
    a.center+a.L*im*(1-y)./(y+1)
end


# algebra
*(c::Number,d::PeriodicLine)=PeriodicLine(isapprox(d.center,0) ? d.center : c*d.center,angle(d)+angle(c))
*(d::PeriodicLine,c::Number)=c*d
for OP in (:+,:-)
    @eval begin
        $OP(c::Number,d::PeriodicLine{a,T}) where {a,T}=PeriodicLine{a,promote_type(eltype(c),T)}($OP(c,d.center),d.L)
        $OP(d::PeriodicLine{a,T},c::Number) where {a,T}=PeriodicLine{a,promote_type(eltype(c),T)}($OP(d.center,c),d.L)
    end
end


arclength(d::Union{Line,PeriodicLine}) = Inf
first(d::Union{Line,PeriodicLine})= -Inf
last(d::Union{Line,PeriodicLine})= Inf
complexlength(d::Union{Line,PeriodicLine})=Inf

## vectorized

for typ in (:Line,:PeriodicLine)
    @eval function convert(::Type{$typ},d::IntervalSets.ClosedInterval)
        a,b=d.left,d.right
        @assert abs(a) == abs(b) == Inf

        if isa(a,Real) && isa(b,Real)
            if a==Inf
                @assert b==-Inf
                $typ(true)
            else
                @assert a==-Inf&&b==Inf
                $typ(false)
            end
        elseif abs(real(a)) < Inf
            @assert real(a)==real(b)
            @assert sign(imag(a))==-sign(imag(b))

            $typ(real(b),angle(b))
        elseif isnan(real(a)) && isnan(real(b))  # hack for -im*Inf
            $typ([imag(a)*im,imag(b)*im])
        elseif isnan(real(a))  # hack for -im*Inf
            $typ([real(b)+imag(a)*im,b])
        elseif isnan(real(b))  # hack for -im*Inf
            $typ([a,real(a)+imag(b)*im])
        elseif abs(imag(a)) < Inf
            @assert imag(a)==imag(b)
            @assert sign(real(a))==-sign(real(b))

            $typ(imag(b),angle(b))
        else
            @assert angle(b) == -angle(a)

            $typ(0.,angle(b))
        end
    end
end
