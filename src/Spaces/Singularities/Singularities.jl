"""
 WeightSpace represents a space that weights another space.
 Overload weight(S,x).
"""
abstract WeightSpace{S,T,DD,d} <: Space{T,DD,d}


domain(S::WeightSpace)=domain(S.space)


points(sp::WeightSpace,n)=points(sp.space,n)


immutable WeightSpacePlan{S,P,T,V}
    space::S
    plan::P
    points::Vector{T}
    weights::Vector{V}
end


for TYP in (:plan_transform,:plan_itransform)
    @eval function $TYP(S::WeightSpace,vals::Vector)
        pts=points(S,length(vals))
        WeightSpacePlan(S,$TYP(S.space,vals),pts,weight(S,pts))
    end
end


transform(sp::WeightSpace,vals::Vector,plan::WeightSpacePlan) =
    transform(sp.space,vals./(sp==plan.space?plan.weights:weight(sp,plan.points)),plan.plan)
itransform(sp::WeightSpace,cfs::Vector,plan::WeightSpacePlan) =
    itransform(sp.space,cfs,plan.plan).*(sp==plan.space?plan.weights:weight(sp,plan.points))



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
include("LogWeight.jl")
include("ExpWeight.jl")
