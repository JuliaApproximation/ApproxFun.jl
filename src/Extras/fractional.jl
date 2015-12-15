####
# Support for Fractional derivatives
#
# defined as
####



export LeftIntegral,LeftDerivative

@calculus_operator LeftIntegral
@calculus_operator LeftDerivative

function rangespace{DD<:Interval}(Q::LeftIntegral{Jacobi{DD},Float64})
    if Q.order==0.5
        S=domainspace(Q)
        @assert S.a==0
        @assert S.b==0
        JacobiWeight(0.5,0.,Jacobi(-0.5,0.5,domain(S)))
    else
        error("Not implemented")
    end
end

bandinds{DD<:Interval}(Q::LeftIntegral{Jacobi{DD},Float64})=(0,0)

jacobi_frac_addentries!(d::Interval,α,μ,A,kr::UnitRange)=
    jacobi_frac_addentries!(sqrt(length(d)/2),α,μ,A,kr)
function jacobi_frac_addentries!(c::Number,α,μ,A,kr::UnitRange)
    γ=c*gamma(α+1)/gamma(α+1+μ)
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
    @assert μ==0.5
    @assert S.a==0
    @assert S.b==0

    # the 1/sqrt(length(d)) gives the constant term
    jacobi_frac_addentries!(domain(S),0.,μ,A,kr)
end


function LeftIntegral{DD}(S::JacobiWeight{Chebyshev{DD}},k)
    # convert to Jacobi
    @assert k==.5

    Q=LeftIntegral(JacobiWeight(S.α,S.β,Jacobi(-.5,.5,domain(S))),k)
    ApproxFun.LeftIntegralWrapper(Q*Conversion(S,domainspace(Q)),k)
end


function rangespace{DD<:Interval}(Q::LeftIntegral{JacobiWeight{Jacobi{DD},DD},Float64})
    μ=Q.order
    S=domainspace(Q)
    J=S.space
    @assert S.β==0
    @assert S.α==J.b
    if isapprox(S.α,-μ)
        Jacobi(J.a-μ,J.b+μ,domain(J))
    else
        JacobiWeight(S.α+μ,0.,Jacobi(J.a-μ,J.b+μ,domain(J)))
    end
end

bandinds{DD<:Interval}(Q::LeftIntegral{JacobiWeight{Jacobi{DD},DD},Float64})=(0,0)

function addentries!{DD<:Interval}(Q::LeftIntegral{JacobiWeight{Jacobi{DD},DD},Float64},A,kr::UnitRange,::Colon)
    μ=Q.order
    S=domainspace(Q)
    J=S.space
    @assert S.β==0
    @assert S.α==J.b

    jacobi_frac_addentries!(domain(S),S.α,μ,A,kr)
end

function choosedomainspace{T<:Float64}(Q::LeftIntegral{UnsetSpace,T},sp::JacobiWeight)
    #TODO: general case
    @assert Q.order==0.5
    @assert isapproxinteger(sp.α-0.5) && sp.α>0 && isapproxinteger(sp.β)
    Legendre(domain(sp))
end

function choosedomainspace{T<:Float64}(Q::LeftIntegral{UnsetSpace,T},sp::PolynomialSpace)
    #TODO: general case
    @assert Q.order==0.5
    JacobiWeight(0.5,0.,Jacobi(0.5,0.5,domain(sp)))
end





function LeftDerivative(S::Space,k)
    i=ceil(Int,k)
    r=i-k
    LeftDerivativeWrapper(i<0?LeftIntegral(S,-k):Derivative(i)*LeftIntegral(S,r),k)
end




## SumSpace

LeftIntegral(S::SumSpace,k)=LeftIntegralWrapper(DiagonalInterlaceOperator(map(s->LeftIntegral(s,k),S.spaces),SumSpace),k)
