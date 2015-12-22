
points(S::Jacobi,n)=fromcanonical(S,gaussjacobi(n,S.a,S.b)[1])

immutable JacobiTransformPlan{DD,T}
    space::Jacobi{DD}
    points::Vector{T}
    weights::Vector{T}
end

plan_transform(S::Jacobi,v::Vector) = JacobiTransformPlan(S,gaussjacobi(length(v),S.a,S.b)...)
function transform(S::Jacobi,vals,plan::JacobiTransformPlan)
#    @assert S==plan.space
    x,w = plan.points, plan.weights
    V=jacobip(0:length(vals)-1,S.a,S.b,x)'
    nrm=(V.^2)*w

    V*(w.*vals)./nrm
end


immutable JacobiITransformPlan{DD,T}
    space::Jacobi{DD}
    points::Vector{T}
end


plan_itransform(S::Jacobi,cfs::Vector) = JacobiITransformPlan(S,points(S,length(cfs)))
function itransform(S::Jacobi,cfs,plan::JacobiITransformPlan)
#    @assert S==plan.space
    jacobip(0:length(cfs)-1,S.a,S.b,tocanonical(S,plan.points))*cfs
end

