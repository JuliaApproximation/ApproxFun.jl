"""
 WeightSpace represents a space that weights another space.
 Overload weight(S,x).
"""
abstract type WeightSpace{S,DD,RR} <: Space{DD,RR} end


domain(S::WeightSpace) = domain(S.space)


points(sp::WeightSpace,n) = points(sp.space,n)


struct WeightSpacePlan{S,P,T,V}
    space::S
    plan::P
    points::Vector{T}
    weights::Vector{V}
end

struct IWeightSpacePlan{S,P,T,V}
    space::S
    plan::P
    points::Vector{T}
    weights::Vector{V}
end

function plan_transform(S::WeightSpace,vals::AbstractVector)
    pts=points(S,length(vals))
    WeightSpacePlan(S,plan_transform(S.space,vals),pts,weight.(S,pts))
end

function plan_itransform(S::WeightSpace,vals::AbstractVector)
    pts=points(S,length(vals))
    IWeightSpacePlan(S,plan_itransform(S.space,vals),pts,weight.(S,pts))
end

*(P::WeightSpacePlan,vals::AbstractVector) = P.plan*(vals./P.weights)
*(P::IWeightSpacePlan,cfs::AbstractVector) = P.weights.*(P.plan*cfs)


# used for ProductFun
transform(sp::WeightSpace,vals::AbstractVector,plan::WeightSpacePlan) =
    transform(sp.space,vals./(sp==plan.space?plan.weights:weight.(sp,plan.points)),plan.plan)
itransform(sp::WeightSpace,cfs::AbstractVector,plan::WeightSpacePlan) =
    itransform(sp.space,cfs,plan.plan).*(sp==plan.space?plan.weights:weight.(sp,plan.points))



function evaluate(f::AbstractVector,S::WeightSpace,x)
    fv=evaluate(f,S.space,x)
    weight(S,x).*fv
end

function evaluate(f::AbstractVector,S::WeightSpace,x...)
    fv=evaluate(f,S.space,x...)
    weight(S,x...).*fv
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
