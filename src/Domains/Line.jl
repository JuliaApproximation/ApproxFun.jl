

export Line, PeriodicLine




## Standard interval


immutable Line{T<:Number} <: IntervalDomain{T}
    centre::T  ##TODO Allow complex
    angle::T
    α::T
    β::T

    #TODO get this inner constructor working again.
    #Line(c,a,α,β)= begin @assert c==a==0.; @assert α<0; @assert β<0; new(c,a,α,β) end
end



Line(c,a)=Line(c,a,-1.,-1.)
Line()=Line(0.,0.)
# true is negative orientation, false is positive orientation
# this is because false==0 and we take angle=0
Line(b::Bool)=b?Line(0.,π):Line()


isambiguous(d::Line)=isnan(d.centre) && isnan(d.angle)
Base.convert{T<:Number}(::Type{Line{T}},::AnyDomain)=Line{T}(NaN,NaN,-1.,-1.)
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


tocanonical(d::Line,x)=line_tocanonical(d.α,d.β,cistyped(-d.angle).*(x-d.centre))
tocanonicalD(d::Line,x)=cistyped(-d.angle).*line_tocanonicalD(d.α,d.β,cistyped(-d.angle).*(x-d.centre))
fromcanonical(d::Line,x)=cistyped(d.angle)*line_fromcanonical(d.α,d.β,x)+d.centre
fromcanonicalD(d::Line,x)=cistyped(d.angle)*line_fromcanonicalD(d.α,d.β,x)
invfromcanonicalD(d::Line,x)=cistyped(-d.angle)*line_invfromcanonicalD(d.α,d.β,x)




Base.length(d::Line) = Inf
Base.first(d::Line)= -Inf
Base.last(d::Line)= Inf


==(d::Line,m::Line) = d.centre == m.centre && d.angle == m.angle && d.β == m.β &&d.α == m.α



# algebra
*(c::Number,d::Line)=Line(isapprox(d.centre,0)?d.centre:c*d.centre,d.angle+angle(c),d.α,d.β)
+(c::Number,d::Line)=Line(c+d.centre,d.angle,d.α,d.β)



## Periodic line

# angle is (false==0) and π (true==1)
# or ranges from (-1,1]
immutable PeriodicLine{angle} <: PeriodicDomain{Float64}
    centre::Float64  ##TODO Allow complex
    L::Float64
    PeriodicLine(c,L)=new(c,L)
    PeriodicLine(c)=new(c,1.)
    PeriodicLine()=new(0.,1.)
end

canonicaldomain(::PeriodicLine)=PeriodicInterval()
PeriodicLine(c,a)=PeriodicLine{a/π}(c,1.)
PeriodicLine()=PeriodicLine{false}(0.,1.)
PeriodicLine(b::Bool)=PeriodicLine{b}()

isambiguous(d::PeriodicLine)=isnan(d.centre) && isnan(d.angle)
Base.convert{T<:Number}(::Type{PeriodicLine{T}},::AnyDomain)=PeriodicLine{T}(NaN,NaN)
Base.convert{IT<:PeriodicLine}(::Type{IT},::AnyDomain)=PeriodicLine(NaN,NaN)

Base.angle{a}(d::PeriodicLine{a})=a*π

tocanonical(d::PeriodicLine{false},x)= 2atan((x-d.centre)/d.L)
fromcanonical(d::PeriodicLine{false},θ)=d.L*tan(θ/2) + d.centre

tocanonical{a}(d::PeriodicLine{a},x)=tocanonical(PeriodicLine{false}(0.,d.L),exp(-π*im*a)*(x-d.centre))
fromcanonical{a}(d::PeriodicLine{a},x)=exp(π*im*a)*fromcanonical(PeriodicLine{false}(0.,d.L),x)+d.centre


function invfromcanonicalD(d::PeriodicLine{false})
    @assert d.centre==0  && d.L==1.0
    a=Fun([1.,0,1],PeriodicInterval())
end

mappoint(a::PeriodicLine{false},b::Circle,x)=b.radius*((a.L*im-(x-a.centre))./(a.L*im+(x-a.centre)))+b.center
function mappoint(b::Circle,a::PeriodicLine{false},x)
    y=(x-b.center)./b.radius
    a.centre+a.L*im*(1-y)./(y+1)
end


# algebra
*(c::Number,d::PeriodicLine)=PeriodicLine(isapprox(d.centre,0)?d.centre:c*d.centre,angle(d)+angle(c))
+{a}(c::Number,d::PeriodicLine{a})=PeriodicLine{a}(c+d.centre,d.L)

## vectorized

for typ in (:Line,:PeriodicLine)
    @eval function ($typ)(d::Vector)
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




