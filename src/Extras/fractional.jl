function rangespace(Q::Integral{Jacobi,Float64})
if Q.order==0.5
    @assert domainspace(Q)==Legendre()
    JacobiWeight(-0.5,0,Chebyshev(domain(Q)))
else
    error("Not implemented")
end
end

function bandinds(Q::Integral{Jacobi,Float64})
   if Q.order==0.5
        (-1,0)
    else
        error("Not implemented")
    end
end

function addentries!(Q::Integral{Jacobi,Float64},A,kr::Range)
    @assert domain(Q)==Interval()
    if Q.order==0.5
        for k=kr
            if k==1
                A[1,1]+=2/sqrt(π)
            else
                A[k,k-1]+=2/((2(k-2)+1)*sqrt(π))
                A[k,k]+=2/((2(k-1)+1)*sqrt(π))
            end
        end
    else
        error("Not implemented")
    end
    A
end

function rangespace(Q::Integral{JacobiWeight{Jacobi},Float64})
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

bandinds(Q::Integral{JacobiWeight{Jacobi},Float64})=(0,0)

function addentries!(Q::Integral{JacobiWeight{Jacobi},Float64},A,kr::UnitRange)
    @assert domain(Q)==Interval()
    μ=Q.order
    S=domainspace(Q)
    J=S.space
    @assert S.β==0
    @assert S.α==J.b


    γ=gamma(S.α)/gamma(S.α+μ)
    for k=1:first(kr)-1
        γ*=(S.α+k-1)/(S.α+μ+k-1)
    end
    for k=kr
       γ*=(S.α+k-1)/(S.α+μ+k-1)
       A[k,k]+=γ
        #should be gamma(S.α+k)/gamma(S.α+μ+k)=
        #(S.α+k-1)/(S.α+μ+k-1)*gamma(S.α+k-1)/gamma(S.α+μ+k-1)
    end
    A
end

function choosedomainspace{T<:Float64}(Q::Integral{UnsetSpace,T},sp::JacobiWeight)
    #TODO: general case
    @assert Q.order==0.5
    @assert sp.α==0.5 && sp.β==0.
    Legendre()⊕JacobiWeight(0.5,0.,Jacobi(0.5,0.5,domain(sp)))
end

function choosedomainspace{T<:Float64}(Q::Integral{UnsetSpace,T},sp::PolynomialSpace)
    #TODO: general case
    @assert Q.order==0.5
    Legendre()⊕JacobiWeight(0.5,0.,Jacobi(0.5,0.5,domain(sp)))
end
