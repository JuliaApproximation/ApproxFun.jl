__precompile__()

module ApproxFun
    using Base, FastGaussQuadrature


import Base: values,getindex,setindex!,*,.*,+,.+,-,.-,==,<,<=,>,
                >=,./,/,.^,^,\,∪,transpose


export pad!,pad,sample,chop!,complexroots,roots,svfft

##Testing
export bisectioninv, clenshaw





include("LinearAlgebra/LinearAlgebra.jl")


include("Fun/Fun.jl")
include("Multivariate/Multivariate.jl")
include("Operators/Operator.jl")

include("Domains/Domains.jl")
include("Spaces/Spaces.jl")




## Further extra features

include("PDE/PDE.jl")
include("Plot/Plot.jl")
include("Extras/Extras.jl")

include("docs.jl")



## precompile
precompile(Chebyshev,tuple())
precompile(Fun,tuple())
precompile(Base.exp,(Fun{Chebyshev{Interval{Float64}},Float64},))
precompile(linsolve,(Vector{Operator{Float64}},Vector{Float64}))
precompile(\,(Vector{Operator{Float64}},Vector{Float64}))
precompile(adaptiveqr,(Vector{Operator{Float64}},Vector{Float64},Float64,Int))
precompile(MutableOperator,(Vector{Operator{Float64}},))
precompile(adaptiveqr!,(MutableOperator,Vector{Float64},Float64,Int))

end #module
