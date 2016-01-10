

export Line, PeriodicLine




## Standard interval



# angle is π*a where a is (false==0) and (true==1)
# or ranges from (-1,1].  We use 1 as 1==true.

immutable Line{angle,T<:Number} <: IntervalDomain{T}
    center::T
    α::Float64
    β::Float64

    #TODO get this inner constructor working again.
    Line(c,α,β)=new(c,α,β)
    Line(c)=new(c,-1.,-1.)
    Line()=new(zero(T),-1.,-1.)
end

typealias RealLine{T} Union{Line{false,T},Line{true,T}}

Base.convert{a}(::Type{Line{a}},c,α,β)=Line{a,typeof(c)}(c,α,β)
Base.convert{a}(::Type{Line{a}},c::Number)=Line{a,typeof(c)}(c)
Base.convert{a}(::Type{Line{a}})=Line{a,Float64}()

Base.angle{a}(d::Line{a})=a*π



Base.reverse(d::Line{true})=Line{false}(d.center,d.β,d.α)
Base.reverse(d::Line{false})=Line{true}(d.center,d.β,d.α)
Base.reverse{a}(d::Line{a})=Line{a-1}(d.center,d.β,d.α)

# ensure the angle is always in (-1,1]
Line(c,a,α,β)=Line{mod(a/π-1,-2)+1,typeof(c)}(c,α,β)
Line(c,a)=Line(c,a,-1.,-1.)
# true is negative orientation, false is positive orientation
# this is because false==0 and we take angle=0
Line(b::Bool)=Line{b}()
Line()=Line(false)

function Line(d::AbstractVector)
    @assert length(d)==2 && isinf(d[1]) && isinf(d[2])

    if d[1]==Inf && d[2] == -Inf
        Line(true)
    elseif d[1]==-Inf && d[2] == Inf
        Line(false)
    else
        error("Not implemented")
    end
end



isambiguous(d::Line)=isnan(d.center)
Base.convert{a,T<:Number}(::Type{Line{a,T}},::AnyDomain)=Line{a,T}(NaN)
Base.convert{IT<:Line}(::Type{IT},::AnyDomain)=Line(NaN,NaN)

## Map interval


##TODO non-1 alpha,beta


function line_tocanonical(α,β,x)


    @assert α==β==-1. || α==β==-.5

    if α==β==-1.
        2x./(1+sqrt(1+4x.^2))
    elseif α==β==-.5
        x./sqrt(1 + x.^2)
    end
end

function line_tocanonicalD(α,β,x)
    @assert α==β==-1. || α==β==-.5

    if α==β==-1.
        2./(1+4x.^2+sqrt(1+4x.^2))
    elseif α==β==-.5
        (1 + x.^2).^(-3/2)
    end
end
function line_fromcanonical(α,β,x)
    #TODO: why is this consistent?
    if α==β==-1.
        x./(1-x.^2)
    else
        x.*(1 + x).^α.*(1 - x).^β
    end
end
function line_fromcanonicalD(α,β,x)
    if α==β==-1.
        (1+x.^2)./(1-x.^2).^2
    else
        (1 - (β-α)x - (β+α+1)x.^2).*(1+x).^(α-1).*(1-x).^(β-1)
    end
end

function line_invfromcanonicalD(α,β,x)
    if α==β==-1.
        (1-x.^2).^2./(1+x.^2)
    else
        1./(1 - (β-α)x - (β+α+1)x.^2).*(1+x).^(1-α).*(1-x).^(1-β)
    end
end


tocanonical(d::Line,x)=line_tocanonical(d.α,d.β,cis(-angle(d)).*(x-d.center))
tocanonical(d::Line{false},x)=line_tocanonical(d.α,d.β,x-d.center)
tocanonical(d::Line{true},x)=line_tocanonical(d.α,d.β,d.center-x)

tocanonicalD(d::Line,x)=cis(-angle(d)).*line_tocanonicalD(d.α,d.β,cis(-angle(d)).*(x-d.center))
tocanonicalD(d::Line{false},x)=line_tocanonicalD(d.α,d.β,x-d.center)
tocanonicalD(d::Line{true},x)=-line_tocanonicalD(d.α,d.β,d.center-x)

fromcanonical(d::Line,x)=cis(angle(d))*line_fromcanonical(d.α,d.β,x)+d.center
fromcanonical(d::Line{false},x)=line_fromcanonical(d.α,d.β,x)+d.center
fromcanonical(d::Line{true},x)=-line_fromcanonical(d.α,d.β,x)+d.center

fromcanonicalD(d::Line,x)=cis(angle(d))*line_fromcanonicalD(d.α,d.β,x)
fromcanonicalD(d::Line{false},x)=line_fromcanonicalD(d.α,d.β,x)
fromcanonicalD(d::Line{true},x)=-line_fromcanonicalD(d.α,d.β,x)

invfromcanonicalD(d::Line,x)=cis(-angle(d))*line_invfromcanonicalD(d.α,d.β,x)
invfromcanonicalD(d::Line{false},x)=line_invfromcanonicalD(d.α,d.β,x)
invfromcanonicalD(d::Line{true},x)=-line_invfromcanonicalD(d.α,d.β,x)






=={a}(d::Line{a},m::Line{a}) = d.center == m.center && d.β == m.β &&d.α == m.α



# algebra
*(c::Real,d::Line{false})=Line{sign(c)>0?false:true}(isapprox(d.center,0)?d.center:c*d.center,d.α,d.β)
*(c::Real,d::Line{true})=Line{sign(c)>0?true:false}(isapprox(d.center,0)?d.center:c*d.center,d.α,d.β)
*(c::Number,d::Line)=Line(isapprox(d.center,0)?d.center:c*d.center,angle(d)+angle(c),d.α,d.β)
*(d::Line,c::Number)=c*d
for OP in (:+,:-)
    @eval begin
        $OP{a}(c::Number,d::Line{a})=Line{a}($OP(c,d.center),d.α,d.β)
        $OP{a}(d::Line{a},c::Number)=Line{a}($OP(d.center,c),d.α,d.β)
    end
end






## Periodic line

# angle is (false==0) and π (true==1)
# or ranges from (-1,1]
immutable PeriodicLine{angle,T} <: PeriodicDomain{Float64}
    center::T
    L::Float64
    PeriodicLine(c,L)=new(c,L)
    PeriodicLine(c)=new(c,1.)
    PeriodicLine()=new(0.,1.)
end

Base.convert{a}(::Type{PeriodicLine{a}},c,L)=PeriodicLine{a,typeof(c)}(c,L)


PeriodicLine(c,a)=PeriodicLine{a/π,eltype(c)}(c,1.)
PeriodicLine()=PeriodicLine{false,Float64}(0.,1.)
PeriodicLine(b::Bool)=PeriodicLine{b,Float64}()

isambiguous(d::PeriodicLine)=isnan(d.center) && isnan(d.angle)
Base.convert{T<:Number,TT}(::Type{PeriodicLine{T,TT}},::AnyDomain)=PeriodicLine{T,TT}(NaN,NaN)
Base.convert{IT<:PeriodicLine}(::Type{IT},::AnyDomain)=PeriodicLine(NaN,NaN)

Base.angle{a}(d::PeriodicLine{a})=a*π

Base.reverse(d::PeriodicLine{true})=PeriodicLine{false}(d.center,d.L)
Base.reverse(d::PeriodicLine{false})=PeriodicLine{true}(d.center,d.L)
Base.reverse{a}(d::PeriodicLine{a})=PeriodicLine{a-1}(d.center,d.L)

tocanonical(d::PeriodicLine{false},x)= 2atan((x-d.center)/d.L)
fromcanonical(d::PeriodicLine{false},θ)=d.L*tan(θ/2) + d.center

tocanonical{a}(d::PeriodicLine{a},x)=tocanonical(PeriodicLine{false,Float64}(0.,d.L),exp(-π*im*a)*(x-d.center))
fromcanonical{a}(d::PeriodicLine{a},x)=exp(π*im*a)*fromcanonical(PeriodicLine{false,Float64}(0.,d.L),x)+d.center


function invfromcanonicalD(d::PeriodicLine{false})
    @assert d.center==0  && d.L==1.0
    a=Fun([1.,0,1],PeriodicInterval())
end

mappoint(a::PeriodicLine{false},b::Circle,x)=b.radius*((a.L*im-(x-a.center))./(a.L*im+(x-a.center)))+b.center
function mappoint(b::Circle,a::PeriodicLine{false},x)
    y=(x-b.center)./b.radius
    a.center+a.L*im*(1-y)./(y+1)
end


# algebra
*(c::Number,d::PeriodicLine)=PeriodicLine(isapprox(d.center,0)?d.center:c*d.center,angle(d)+angle(c))
*(d::PeriodicLine,c::Number)=c*d
for OP in (:+,:-)
    @eval begin
        $OP{a,T}(c::Number,d::PeriodicLine{a,T})=PeriodicLine{a,promote_type(eltype(c),T)}($OP(c,d.center),d.L)
        $OP{a,T}(d::PeriodicLine{a,T},c::Number)=PeriodicLine{a,promote_type(eltype(c),T)}($OP(d.center,c),d.L)
    end
end


Base.length(d::Union{Line,PeriodicLine}) = Inf
Base.first(d::Union{Line,PeriodicLine})= -Inf
Base.last(d::Union{Line,PeriodicLine})= Inf
complexlength(d::Union{Line,PeriodicLine})=Inf

## vectorized

for typ in (:Line,:PeriodicLine)
    @eval function ($typ)(d::AbstractVector)
        @assert length(d) ==2
        @assert abs(d[1]) == abs(d[2]) == Inf

        if isa(d[1],Real) && isa(d[2],Real)
            if d[1]==Inf
                @assert d[2]==-Inf
                $typ(true)
            else
                @assert d[1]==-Inf&&d[2]==Inf
                $typ(false)
            end
        elseif abs(real(d[1])) < Inf
            @assert real(d[1])==real(d[2])
            @assert sign(imag(d[1]))==-sign(imag(d[2]))

            $typ(real(d[2]),angle(d[2]))
        elseif isnan(real(d[1])) && isnan(real(d[2]))  # hack for -im*Inf
            $typ([imag(d[1])*im,imag(d[2])*im])
        elseif isnan(real(d[1]))  # hack for -im*Inf
            $typ([real(d[2])+imag(d[1])*im,d[2]])
        elseif isnan(real(d[2]))  # hack for -im*Inf
            $typ([d[1],real(d[1])+imag(d[2])*im])
        elseif abs(imag(d[1])) < Inf
            @assert imag(d[1])==imag(d[2])
            @assert sign(real(d[1]))==-sign(real(d[2]))

            $typ(imag(d[2]),angle(d[2]))
        else
            @assert angle(d[2]) == -angle(d[1])

            $typ(0.,angle(d[2]))
        end
    end
end
