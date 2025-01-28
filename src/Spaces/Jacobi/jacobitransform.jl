
points(S::Jacobi,n)=fromcanonical(S,gaussjacobi(n,S.a,S.b)[1])

immutable JacobiTransformPlan{DD,T,TT}
    space::Jacobi{TT,DD}
    points::Vector{T}
    weights::Vector{T}
end

plan_transform(S::Jacobi,v::Vector) = JacobiTransformPlan(S,gaussjacobi(length(v),S.a,S.b)...)
function transform(S::Jacobi,vals,plan::JacobiTransformPlan)
#    @assert S==plan.space
    x,w = plan.points, plan.weights
    V=jacobip(0:length(vals)-1,S.a,S.b,x)'
    w2=w.*(1-x).^(S.a-plan.space.a).*(1+x).^(S.b-plan.space.b)   # need to weight if plan is different
    nrm=(V.^2)*w2

    V*(w2.*vals)./nrm
end




immutable JacobiITransformPlan{DD,T}
    space::Jacobi{Float64,DD}
    points::Vector{T}
end


plan_itransform(S::Jacobi,cfs::Vector) = JacobiITransformPlan(S,points(S,length(cfs)))
function itransform(S::Jacobi,cfs,plan::JacobiITransformPlan)
#    @assert S==plan.space
    jacobip(0:length(cfs)-1,S.a,S.b,tocanonical(S,plan.points))*cfs
end

function coefficients(f::Vector,a::Jacobi,b::Chebyshev)
    if domain(a) == domain(b) && (!isapproxinteger(a.a-0.5) || !isapproxinteger(a.b-0.5))
        cjt(f,a.a,a.b)
    else
        defaultcoefficients(f,a,b)
    end
end
function coefficients(f::Vector,a::Chebyshev,b::Jacobi)
    if domain(a) == domain(b) && (!isapproxinteger(b.a-0.5) || !isapproxinteger(b.b-0.5))
        icjt(f,b.a,b.b)
    else
        defaultcoefficients(f,a,b)
    end
end

function coefficients(f::Vector,a::Jacobi,b::Jacobi)
    if domain(a) == domain(b) && (!isapproxinteger(a.a-b.a) || !isapproxinteger(a.b-b.b))
        jjt(f,a.a,a.b,b.a,b.b)
    else
        defaultcoefficients(f,a,b)
    end
end
