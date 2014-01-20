

export Ray



## Standard interval


type Ray <: IntervalDomain
    centre::Float64  ##TODO Allow complex
    angle::Float64
    
    Ray(c,a)=(@assert c==a==0.; new(c,a))
end

Ray()=Ray(0.,0.)

## Map interval



tocanonical(d::Ray,x)=(x-1.)./(1.+x)
tocanonicalD(d::Ray,x)=2./((1.+x).^2)
fromcanonical(d::Ray,x)=(1.+x)./(1.-x)
fromcanonicalD(d::Ray,x)=2./((x-1.).^2)



Base.length(d::Ray) = Inf



==(d::Ray,m::Ray) = d.centre == m.centre && d.angle == m.angle



