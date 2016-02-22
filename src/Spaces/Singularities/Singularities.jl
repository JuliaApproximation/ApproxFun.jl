"""
 WeightSpace represents a space that weights another space.
 Overload weight(S,x).
"""
abstract WeightSpace{T,DD,d} <: Space{T,DD,d}


domain(S::WeightSpace)=domain(S.space)


points(sp::WeightSpace,n)=points(sp.space,n)
plan_transform(S::WeightSpace,vals::Vector)=1./weight(S,points(S,length(vals))),plan_transform(S.space,vals)
plan_itransform(S::WeightSpace,vals::Vector)=weight(S,points(S,length(vals))),plan_itransform(S.space,vals)


transform(sp::WeightSpace,vals::Vector,plan)=transform(sp.space,vals.*plan[1],plan[2])
itransform(sp::WeightSpace,cfs::Vector,plan)=itransform(sp.space,cfs,plan[2]).*plan[1]



function evaluate(f::AbstractVector,S::WeightSpace,x...)
    tol=1.0E-14
    fv=Fun(f,S.space)(x...)
    if isa(fv,Number)&&abs(fv)<tol
        #TODO: Why this special case??
        zero(eltype(fv))
    else
        weight(S,x...).*fv
    end
end




include("divide_singularity.jl")
include("JacobiWeight.jl")
include("JacobiWeightOperators.jl")
include("JacobiWeightChebyshev.jl")
include("HeavisideSpace.jl")
include("DiracSpace.jl")
include("PointSpace.jl")
include("LogWeight.jl")
