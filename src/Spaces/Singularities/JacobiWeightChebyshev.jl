

function coefficients{C<:Chebyshev,DD<:IntervalDomain}(f::Vector,sp1::JacobiWeight{C,DD},sp2::JacobiWeight{C,DD})
    α,β=sp1.α,sp1.β
    c,d=sp2.α,sp2.β
    tol=10E-5  #TODO:tol choice
    @assert isapproxinteger(α-c) && isapproxinteger(β-d)
    if isapprox(c,α) && isapprox(d,β)
        f
    elseif c>α && d>β
        g=Fun(f,sp1.space)
        @assert abs(first(g))<tol && abs(last(g))<tol
        coefficients(divide_singularity(g).coefficients,JacobiWeight(α+1,β+1,sp1.space),sp2)
    elseif c>α
        g=Fun(f,sp1.space)
        @assert abs(first(g))<tol
        coefficients(divide_singularity(false,g).coefficients,JacobiWeight(α+1,β,sp1.space),sp2)
    elseif d>β
        g=Fun(f,sp1.space)
        @assert abs(last(g))<tol
        coefficients(divide_singularity(true,g).coefficients,JacobiWeight(α,β+1,sp1.space),sp2)
    elseif c<α && d<β
        x=Fun()
        coefficients(((1-x^2)*Fun(f,Chebyshev())).coefficients,JacobiWeight(α-1,β-1,sp1.space),sp2)
    elseif d<β
        x=Fun()
        coefficients(((1-x)*Fun(f,Chebyshev())).coefficients,JacobiWeight(α,β-1,sp1.space),sp2)
    elseif c<α
        x=Fun()
        coefficients(((1+x)*Fun(f,Chebyshev())).coefficients,JacobiWeight(α-1,β,sp1.space),sp2)
    else
        error("Something has gone wrong")
    end
end
