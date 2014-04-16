

export Ray



## Standard interval


type Ray <: IntervalDomain
    centre::Float64  ##TODO Allow complex
    angle::Float64
    
    Ray(c,a)=(@assert a==0.; new(c,a))
end

Ray()=Ray(0.,0.)

## Map interval



ray_tocanonical(x)=(x.-1.)./(1.+x)
ray_tocanonicalD(x)=2./((1.+x).^2)
ray_fromcanonical(x)=(1.+x)./(1.-x)
ray_fromcanonicalD(x)=2./((x.-1.).^2)

tocanonical(d::Ray,x)=ray_tocanonical(x.-d.centre)
tocanonicalD(d::Ray,x)=ray_tocanonicalD(x.-d.centre)
fromcanonical(d::Ray,x)=ray_fromcanonical(x).+d.centre
fromcanonicalD(d::Ray,x)=ray_fromcanonicalD(x)




Base.length(d::Ray) = Inf



==(d::Ray,m::Ray) = d.centre == m.centre && d.angle == m.angle





