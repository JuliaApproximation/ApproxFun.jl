"""
 WeightSpace represents a space that weights another space.
 Overload weight(S,x).
"""
abstract WeightSpace{S,T,DD,d} <: Space{T,DD,d}


domain(S::WeightSpace) = domain(S.space)


points(sp::WeightSpace,n) = points(sp.space,n)


immutable WeightSpacePlan{S,P,T,V}
    space::S
    plan::P
    points::Vector{T}
    weights::Vector{V}
end

immutable IWeightSpacePlan{S,P,T,V}
    space::S
    plan::P
    points::Vector{T}
    weights::Vector{V}
end

function plan_transform(S::WeightSpace,vals::Vector)
    pts=points(S,length(vals))
    WeightSpacePlan(S,plan_transform(S.space,vals),pts,weight(S,pts))
end

function plan_itransform(S::WeightSpace,vals::Vector)
    pts=points(S,length(vals))
    IWeightSpacePlan(S,plan_itransform(S.space,vals),pts,weight(S,pts))
end

*(P::WeightSpacePlan,vals::Vector) = P.plan*(vals./P.weights)
*(P::IWeightSpacePlan,cfs::Vector) = P.weights.*(P.plan*cfs)


# transform(sp::WeightSpace,vals::Vector,plan::WeightSpacePlan) =
#     transform(sp.space,vals./(sp==plan.space?plan.weights:weight(sp,plan.points)),plan.plan)
# itransform(sp::WeightSpace,cfs::Vector,plan::WeightSpacePlan) =
#     itransform(sp.space,cfs,plan.plan).*(sp==plan.space?plan.weights:weight(sp,plan.points))



function evaluate(f::AbstractVector,S::WeightSpace,x...)
    tol=1.0E-14
    fv=Fun(S.space,f)(x...)
    if isa(fv,Number)&&abs(fv)<tol
        #TODO: Why this special case??
        zero(eltype(fv))
    else
        weight(S,x...).*fv
    end
end

# recurrence is inhereted
for FUNC in (:recα,:recβ,:recγ)
    @eval $FUNC(T,ws::WeightSpace,k) = $FUNC(T,ws.space,k)
end


include("divide_singularity.jl")
include("JacobiWeight.jl")
include("JacobiWeightOperators.jl")
include("JacobiWeightChebyshev.jl")
include("HeavisideSpace.jl")
include("DiracSpace.jl")
include("LogWeight.jl")
include("ExpWeight.jl")
