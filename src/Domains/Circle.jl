

export Circle


##  Circle


type Circle{T<:Number} <: PeriodicDomain
	center{T}
	radius::Float64
end

Circle(r)=Circle(0.,r)
Circle()=Circle(1.)


function tocanonical(d::Circle,ζ)
    v=(ζ-d.center)/d.radius
    v==-1.0 ? -1.π : angle(v)
end

tocanonicalD(d::Circle,ζ)=-1.im./(ζ-d.center)  #TODO: Check formula
fromcanonical(d::Circle,θ)=d.radius*exp(1.im*θ) + d.center
fromcanonicalD(d::Circle,θ)=d.radius*1.im*exp(1.im*θ)



Base.length(d::Circle) = 2π*d.radius



==(d::Circle,m::Circle) = d.center == m.center && d.radius == m.radius




## Diff and integration


function Base.diff{T<:Number}(f::FFun{T,Circle}) 
        ##TODO: general radii
        @assert f.domain.radius == 1.
        @assert f.domain.center == 0
        cfs = f.coefficients
        # Now shift everything by one
        FFun(ShiftVector(
                        [([firstindex(cfs):-1].*cfs[firstindex(cfs):-1]),0],
                        [1:lastindex(cfs)].*cfs[1:lastindex(cfs)]
                        ),
            f.domain)    
end



function integrate{T<:Number}(f::FFun{T,Circle}) 
    tol = 10eps()
    @assert abs(f.coefficients[-1]) < tol        
    ##TODO: general radii        
    @assert f.domain.radius == 1.
    @assert f.domain.center == 0        
    
    cfs = f.coefficients
    # Now shift everything by one
    FFun(ShiftVector(
                    [cfs[firstindex(cfs):-1]./[firstindex(cfs):-1]],
                    [0,(cfs[0:lastindex(cfs)]./[1:lastindex(cfs)+1])]
                    ),
        f.domain)
end

