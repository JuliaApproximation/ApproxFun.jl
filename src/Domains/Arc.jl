export Arc

immutable Arc{T<:Number,V<:Real} <: IntervalDomain{Complex{Float64}}
    center::T
    radius::V
    angles::(V,V)
end


Arc(c,r,t)=Arc{typeof(c),typeof(r)}(c,r,t)
Arc(c,r,t0,t1)=Arc(c,r,(t0,t1))

intervaltoarc(z0,r,t0,t1,x)=(-exp(im/2*(t0+2t1))*r*(1+x) - 
     exp(im*t0/2)*(1+x)*z0 + 
     exp(im*t1/2)*(x-1)*(exp(im*t0)*r+z0))./(exp(im*t1/2)*(x-1)-exp(im*t0/2)*(1+x))
arctointerval(z0,r,t0,t1,z)=-((exp(im*t0/2)+exp(im*t1/2))*(exp(im*t0/2+im*t1/2)*r-z+ z0))/
((-exp(im*t0/2) + exp(im*t1/2))*(exp(im*t0/2+im*t1/2)*r+z-z0))

fromcanonical(a::Arc,x)=intervaltoarc(a.center,a.radius,a.angles...,x)
tocanonical(a::Arc,x)=arctointerval(a.center,a.radius,a.angles...,x)
#Base.first(a::Arc)=fromcanonical(a,-1.0)


