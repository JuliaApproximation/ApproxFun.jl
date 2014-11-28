

export Line, PeriodicLine




## Standard interval


immutable Line <: IntervalDomain
    centre::Float64  ##TODO Allow complex
    angle::Float64
    α::Float64
    β::Float64    
    
    Line(c,a,α,β)=(@assert c==a==0.; @assert α<0; @assert β<0; new(c,a,α,β))
end



Line(c,a)=Line(c,a,-1.,-1.)
Line()=Line(0.,0.)


## Map interval


##TODO non-1 alpha,beta

function tocanonical(d::Line,x::Number)
    @assert d.α==d.β==-1. || d.α==d.β==-.5
    
    if d.α==d.β==-1.
        abs(x) < 10eps()?0.:.5(sqrt(1+4x.^2) - 1)./x
    elseif d.α==d.β==-.5
        x./sqrt(1 + x.^2)
    end
end

function tocanonical(d::Line,v::Vector)
    @assert d.α==d.β==-1. || d.α==d.β==-.5
    
    if d.α==d.β==-1.
        map(x->abs(x) < 10eps()?0.:.5(sqrt(1+4x.^2) - 1)./x,v)
    elseif d.α==d.β==-.5
        v./sqrt(1 + v.^2)
    end
end
function tocanonicalD(d::Line,v::Vector)
    @assert d.α==d.β==-1. || d.α==d.β==-.5
    
    if d.α==d.β==-1.
        map(x->(abs(x) < 10eps()?
            1.:
            -.5((1 + 4x.^2).^(-.5) - 1)./x.^2
        ),v)
    elseif d.α==d.β==-.5
        (1 + v.^2).^(-3/2)
    end    
end
fromcanonical(d::Line,x)=x.*(1 + x).^d.α.*(1 - x).^d.β
fromcanonicalD(d::Line,x)=(1 - (d.β-d.α)x - (d.β+d.α+1)x.^2).*(1+x).^(d.α-1).*(1-x).^(d.β-1)



Base.length(d::Line) = Inf
Base.first(d::Line)= -Inf
Base.last(d::Line)= Inf


==(d::Line,m::Line) = d.centre == m.centre && d.angle == m.angle && d.β == m.β &&d.α == m.α



##multiplybyx




# function multiplybyx{T<:Number,D<:LineSpace}(f::Fun{T,D})
#     ct=Fun(x->x.*cot(π*x/2),28)
#     x=Fun(identity)
#     u=SingFun(ct./(1-x.^2),1.,1.)
#     Fun((x.*Fun(f)./u).fun./(1-x.^2),domain(f))
# end







## Periodic line

immutable PeriodicLine <: PeriodicDomain
    centre::Float64  ##TODO Allow complex
    angle::Float64
    L::Float64
end

PeriodicLine(c,a)=PeriodicLine(c,a,1.)
PeriodicLine()=PeriodicLine(0.,0.)



tocircle(d::PeriodicLine,x)=(d.L*im - exp(-im*d.angle)*(x-d.centre))./(d.L*im + exp(-im*d.angle)*(x-d.centre))
fromcircle(d::PeriodicLine,ζ)=exp(im*d.angle)*1.im*d.L*(ζ - 1)./(ζ + 1)+d.centre

tocanonical(d::PeriodicLine,x)=tocanonical(Circle(),tocircle(d,x))
fromcanonical(d::PeriodicLine,θ)=fromcircle(d,fromcanonical(Circle(),θ))





## vectorized

for typ in (:Line,:PeriodicLine)
    @eval function ($typ)(d::Vector)
        @assert length(d) ==2
        @assert abs(d[1]) == abs(d[2]) == Inf
        
        if abs(real(d[1])) < Inf 
            @assert real(d[1])==real(d[2])
            @assert sign(imag(d[1]))==-sign(imag(d[2]))
            
            $typ(real(d[2]),angle(d[2]))
            
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







