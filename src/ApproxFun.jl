module ApproxFun
    using Base

export AbstractFun, Fun,FFun,IFun,Interval,evaluate,values,points,chebyshevtransform
export pad!,pad,sample,chop!,complexroots,roots,svfft
export multiplybyx,IntervalDomain,fasttimes

##Testing
export bisectioninv, clenshaw
export coefficients, integrate

export domain

import Base.values

##TODO Incorporate type
abstract AbstractFun

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

FFun(x,d::PeriodicDomain)=IFun(x,LaurentSpace(d))
FFun(x,d::PeriodicDomain,n...)=IFun(x,LaurentSpace(d),n...)
FFun(f,n::Integer)=IFun(f,LaurentSpace(PeriodicInterval()),n)
FFun(f)=IFun(f,LaurentSpace(PeriodicInterval()))

Fun(x)=IFun(x)
Fun(x,d::IntervalDomain)=IFun(x,d)
Fun(x,d::PeriodicDomain)=FFun(x,d)
Fun(x,d)=IFun(x,d)
Fun(x,d::IntervalDomain,n::Integer)=IFun(x,d,n)
Fun(x,d::PeriodicDomain,n::Integer)=FFun(x,d,n)
Fun(x,d::Vector,n::Integer)=IFun(x,d,n)


## General routines

domain(f::IFun)=domain(f.space)
#domain(f::FFun)=f.domain
domain(::Number)=Any
domain{T<:AbstractFun}(v::Vector{T})=map(domain,v)



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


