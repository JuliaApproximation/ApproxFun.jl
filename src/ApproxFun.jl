module ApproxFun
    using Base

export Fun,FFun,Interval,evaluate,values,points,chebyshevtransform
export pad!,pad,sample,chop!,complexroots,roots,svfft
export multiplybyx,IntervalDomain,fasttimes

##Testing
export bisectioninv, clenshaw
export coefficients, integrate

export domain

import Base.values



include("LinearAlgebra/LinearAlgebra.jl")

include("Domains/Domain.jl")
include("Spaces/FunctionSpace.jl")


##Chebyshev Routines
include("IFun/IFun.jl")


# Canonical domains
include("Domains/Interval.jl")
include("Domains/PeriodicInterval.jl")
include("Domains/Ray.jl")
include("Domains/Circle.jl")
include("Spaces/Spaces.jl")

include("Operators/Operator.jl")


include("Multivariate/Multivariate.jl")




##Following routine decides
# whether input is IFun or FFun

FFun(x,d::PeriodicDomain)=Fun(x,LaurentSpace(d))
FFun(x,d::PeriodicDomain,n...)=Fun(x,LaurentSpace(d),n...)
FFun(f,n::Integer)=Fun(f,LaurentSpace(PeriodicInterval()),n)
FFun(f)=Fun(f,LaurentSpace(PeriodicInterval()))




## General routines

domain(f::Fun)=domain(f.space)
#domain(f::FFun)=f.domain
domain(::Number)=Any
domain{T<:Fun}(v::Vector{T})=map(domain,v)



## Other domains


## Further extra features

include("Singularities/Singularities.jl")

include("PDE/PDE.jl")


include("Plot/Plot.jl")


include("sample.jl")





function Interval{T<:Number}(d::Vector{T})
    @assert length(d) >1

    if length(d) == 2    
        if abs(d[1]) == Inf && abs(d[2]) == Inf
            Line(d)
        elseif abs(d[2]) == Inf || abs(d[1]) == Inf
            Ray(d)
        else
            Interval(d[1],d[2])
        end
    else
        [Interval(d[1:2]),Interval(d[2:end])]
    end
end



function PeriodicInterval{T<:Number}(d::Vector{T})
    @assert length(d) == 2
    
    if abs(d[1]) ==Inf
        PeriodicLine(d)
    else
        PeriodicInterval(d[1],d[2])
    end
end






end #module


