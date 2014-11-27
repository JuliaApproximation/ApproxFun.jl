

export Ray



## Standard interval


immutable Ray <: IntervalDomain
    centre::Float64  
    angle::Float64
    orientation::Bool
end

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




## Map interval



ray_tocanonical(x)=(x==Inf)?1.:(x-1.)./(1+x)
ray_tocanonicalD(x)=(x==Inf)?0.:2*(1./(1+x)).^2
ray_fromcanonical(x)=(1+x)./(1-x)
ray_fromcanonicalD(x)=2*(1./(x-1.)).^2

ray_tocanonical(o,x)=(o?1:-1)*ray_tocanonical(x)
ray_tocanonicalD(o,x)=(o?1:-1)*ray_tocanonical(x)
ray_fromcanonical(o,x)=ray_fromcanonical((o?1:-1)*x)
ray_fromcanonicalD(o,x)=(o?1:-1)*ray_fromcanonicalD((o?1:-1)*x)

function cistyped(a::Real)
    if a==0
        1.
    elseif a==π||a==-π
        -1.
    else
        cis(a)
    end
end

tocanonical(d::Ray,x)=ray_tocanonical(d.orientation,cistyped(-d.angle).*(x-d.centre))
tocanonicalD(d::Ray,x)=cistyped(-d.angle).*ray_tocanonicalD(d.orientation,cistyped(-d.angle).*(x-d.centre))
fromcanonical(d::Ray,x)=cistyped(d.angle)*ray_fromcanonical(d.orientation,x)+d.centre
fromcanonicalD(d::Ray,x)=cistyped(d.angle)*ray_fromcanonicalD(d.orientation,x)




Base.length(d::Ray) = Inf



==(d::Ray,m::Ray) = d.centre == m.centre && d.angle == m.angle





