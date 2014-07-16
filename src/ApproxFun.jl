module ApproxFun
    using Base, PyPlot

export Fun,IFun,FFun,Interval,evaluate,values,points,chebyshevtransform
export pad!,pad,sample,chop!,complexroots,roots,svfft
export multiplybyx,IntervalDomain,fasttimes

##Testing
export bisectioninv, clenshaw
export coefficients, integrate

export domain

import PyPlot.plot
export plot


abstract AbstractFun

include("LinearAlgebra/LinearAlgebra.jl")

include("Domains/Domain.jl")



##Chebyshev Routines
include("IFun/IFun.jl")


##Fourier Routines
include("FFun/FFun.jl")

# Canonical domains
include("Domains/Interval.jl")
include("Domains/PeriodicInterval.jl")


include("Operators/Operator.jl")


include("Multivariate/VectorFun.jl")
include("Multivariate/Fun2D.jl")

include("sample.jl")


##Following routine decides
# whether input is IFun or FFun
Fun(x)=IFun(x)
Fun(x,d::IntervalDomain)=IFun(x,d)
Fun(x,d::PeriodicDomain)=FFun(x,d)
Fun(x,d)=IFun(x,d)
Fun(x,d::IntervalDomain,n::Integer)=IFun(x,d,n)
Fun(x,d::PeriodicDomain,n::Integer)=FFun(x,d,n)
Fun(x,d::Vector,n::Integer)=IFun(x,d,n)


## General routines

domain(f::AbstractFun)=f.domain
domain(::Number)=Any
domain{T<:AbstractFun}(v::Vector{T})=map(domain,v)



## Other domains

include("Domains/Ray.jl")
include("Domains/Line.jl")
include("Domains/Circle.jl")

## Further extra features

include("Singularities/SingFun.jl")

include("PDE/pdesolve.jl")

include("Plotting/Plot.jl")



function Interval{T<:Number}(d::Vector{T})
    @assert length(d) >1

    if length(d) == 2    
        if d[1] ==-Inf && d[2] == Inf
            Line()
        elseif d[2] == Inf
            Ray(0,0)
        elseif d[1] == -Inf
            Ray(0,π)
        else
            Interval(d[1],d[2])
        end
    else
        [Interval(d[1],d[2]),Interval(d[2:end])]
    end
end



function PeriodicInterval(d::Vector)
    @assert length(d) == 2
    
    if abs(d[1]) ==Inf
        PeriodicLine(d)
    else
        PeriodicInterval(d[1],d[2])
    end
end






end #module


