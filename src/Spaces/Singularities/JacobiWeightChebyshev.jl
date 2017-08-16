

function coefficients(f::AbstractVector,sp1::JacobiWeight{C,DD},sp2::JacobiWeight{C,DD}) where {C<:Chebyshev,DD<:IntervalDomain}
    β,α=sp1.β,sp1.α
    c,d=sp2.β,sp2.α
    tol=10E-5  #TODO:tol choice
    @assert isapproxinteger(β-c) && isapproxinteger(α-d)
    if isapprox(c,β) && isapprox(d,α)
        f
    elseif c>β && d>α
        g=Fun(sp1.space,f)
        @assert abs(first(g))<tol && abs(last(g))<tol
        coefficients(divide_singularity(g).coefficients,JacobiWeight(β+1,α+1,sp1.space),sp2)
    elseif c>β
        g=Fun(sp1.space,f)
        @assert abs(first(g))<tol
        coefficients(divide_singularity(false,g).coefficients,JacobiWeight(β+1,α,sp1.space),sp2)
    elseif d>α
        g=Fun(sp1.space,f)
        @assert abs(last(g))<tol
        coefficients(divide_singularity(true,g).coefficients,JacobiWeight(β,α+1,sp1.space),sp2)
    elseif c<β && d<α
        x=Fun()
        coefficients(((1-x^2)*Fun(Chebyshev(),f)).coefficients,JacobiWeight(β-1,α-1,sp1.space),sp2)
    elseif d<α
        x=Fun()
        coefficients(((1-x)*Fun(Chebyshev(),f)).coefficients,JacobiWeight(β,α-1,sp1.space),sp2)
    elseif c<β
        x=Fun()
        coefficients(((1+x)*Fun(Chebyshev(),f)).coefficients,JacobiWeight(β-1,α,sp1.space),sp2)
    else
        error("Something has gone wrong")
    end
end
