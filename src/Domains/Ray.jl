

export Ray



## Standard interval


type Ray <: IntervalDomain
    centre::Float64  ##TODO Allow complex
    angle::Float64
end

Ray()=Ray(0.,0.)

## Map interval



ray_tocanonical(x)=(x.-1.)./(1.+x)
ray_tocanonicalD(x)=2./((1.+x).^2)
ray_fromcanonical(x)=(1.+x)./(1.-x)
ray_fromcanonicalD(x)=2./((x.-1.).^2)

function cistyped(a::Real)
    if a==0
        1.
    elseif a==π
        -1.
    else
        cis(a)
    end
end

tocanonical(d::Ray,x)=ray_tocanonical(cistyped(-d.angle).*(x.-d.centre))
tocanonicalD(d::Ray,x)=cistyped(-d.angle).*ray_tocanonicalD(cistyped(-d.angle).*(x.-d.centre))
fromcanonical(d::Ray,x)=cistyped(d.angle)*ray_fromcanonical(x).+d.centre
fromcanonicalD(d::Ray,x)=cistyped(d.angle)*ray_fromcanonicalD(x)




Base.length(d::Ray) = Inf



==(d::Ray,m::Ray) = d.centre == m.centre && d.angle == m.angle





