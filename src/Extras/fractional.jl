@calculus_operator LeftIntegral

function rangespace{DD<:Interval}(Q::LeftIntegral{Jacobi{DD},Float64})
    if Q.order==0.5
        @assert domainspace(Q)==Legendre()
        JacobiWeight(0.5,0.,Jacobi(-0.5,0.5))
    else
        error("Not implemented")
    end
end

bandinds{DD<:Interval}(Q::LeftIntegral{Jacobi{DD},Float64})=(0,0)

function jacobi_frac_addentries!(α,μ,A,kr::UnitRange)
    γ=gamma(α+1)/gamma(α+1+μ)
    for k=1:first(kr)-1
        γ*=(α+k)/(α+μ+k)
    end
    for k=kr
       A[k,k]+=γ
       γ*=(α+k)/(α+μ+k)
        #should be gamma(S.α+k)/gamma(S.α+μ+k)=
        #(S.α+k-1)/(S.α+μ+k-1)*gamma(S.α+k-1)/gamma(S.α+μ+k-1)
    end
    A
end


function addentries!{DD<:Interval}(Q::LeftIntegral{Jacobi{DD},Float64},A,kr::UnitRange,::Colon)
    μ=Q.order
    S=domainspace(Q)
    @assert S==Legendre()

    jacobi_frac_addentries!(0.,μ,A,kr)
end

function rangespace{DD<:Interval}(Q::LeftIntegral{JacobiWeight{Jacobi{DD},DD},Float64})
    μ=Q.order
    S=domainspace(Q)
    J=S.space
    @assert S.β==0
    @assert S.α==J.b
    if isapprox(S.α,-μ)
        Jacobi(J.a-μ,J.b+μ)
    else
        JacobiWeight(S.α+μ,0.,Jacobi(J.a-μ,J.b+μ))
    end
end

bandinds{DD<:Interval}(Q::LeftIntegral{JacobiWeight{Jacobi{DD},DD},Float64})=(0,0)

function addentries!{DD<:Interval}(Q::LeftIntegral{JacobiWeight{Jacobi{DD},DD},Float64},A,kr::UnitRange,::Colon)
    @assert domain(Q)==Interval()
    μ=Q.order
    S=domainspace(Q)
    J=S.space
    @assert S.β==0
    @assert S.α==J.b

    jacobi_frac_addentries!(S.α,μ,A,kr)
end

function choosedomainspace{T<:Float64}(Q::LeftIntegral{UnsetSpace,T},sp::JacobiWeight)
    #TODO: general case
    @assert Q.order==0.5
    @assert sp.α==0.5 && sp.β==0.
    Legendre()⊕JacobiWeight(0.5,0.,Jacobi(0.5,0.5,domain(sp)))
end

function choosedomainspace{T<:Float64}(Q::LeftIntegral{UnsetSpace,T},sp::PolynomialSpace)
    #TODO: general case
    @assert Q.order==0.5
    Legendre()⊕JacobiWeight(0.5,0.,Jacobi(0.5,0.5,domain(sp)))
end
