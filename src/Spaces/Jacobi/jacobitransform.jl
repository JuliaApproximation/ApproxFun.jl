
points(S::Jacobi, n) = points(Chebyshev(domain(S)), n)

struct JacobiTransformPlan{CPLAN,CJT}
    chebplan::CPLAN
    cjtplan::CJT
end

plan_transform(S::Jacobi, v::AbstractVector) =
    JacobiTransformPlan(plan_transform(Chebyshev(), v), plan_icjt(v, S.a, S.b))
*(P::JacobiTransformPlan, vals) = P.cjtplan*(P.chebplan*vals)


struct JacobiITransformPlan{CPLAN,CJT}
    ichebplan::CPLAN
    icjtplan::CJT
end


plan_itransform(S::Jacobi,v::AbstractVector) =
    JacobiITransformPlan(plan_itransform(Chebyshev(), v), plan_cjt(v, S.a, S.b))
*(P::JacobiITransformPlan,cfs) = P.ichebplan*(P.icjtplan*cfs)


function coefficients(f::AbstractVector,a::Jacobi,b::Chebyshev)
    if domain(a) == domain(b) && (!isapproxinteger(a.a-0.5) || !isapproxinteger(a.b-0.5))
        cjt(f,a.a,a.b)
    else
        defaultcoefficients(f,a,b)
    end
end
function coefficients(f::AbstractVector,a::Chebyshev,b::Jacobi)
    isempty(f) && return f
    if domain(a) == domain(b) && (!isapproxinteger(b.a-0.5) || !isapproxinteger(b.b-0.5))
        icjt(f,b.a,b.b)
    else
        defaultcoefficients(f,a,b)
    end
end

function coefficients(f::AbstractVector,a::Jacobi,b::Jacobi)
    if domain(a) == domain(b) && (!isapproxinteger(a.a-b.a) || !isapproxinteger(a.b-b.b))
        jjt(f,a.a,a.b,b.a,b.b)
    else
        defaultcoefficients(f,a,b)
    end
end
