####
# Support for Fractional derivatives
#
# defined as
####



export LeftIntegral,LeftDerivative, RightDerivative, RightIntegral

@calculus_operator LeftIntegral
@calculus_operator LeftDerivative

@calculus_operator RightIntegral
@calculus_operator RightDerivative

# DLMF18.17.9 with μ=0.5 and α=β=0

function LeftIntegral(S::Jacobi,k)
    if S.b==0
        ConcreteLeftIntegral(S,k)
    else
        J=Jacobi(S.a,0.,domain(S))
        LeftIntegralWrapper(LeftIntegral(J,k)*Conversion(S,J),k)
    end
end

function LeftIntegral(S::Ultraspherical,k)
    if order(S) == 1/2
        LeftIntegralWrapper(SpaceOperator(LeftIntegral(Jacobi(S),k),S,S),0.5)
    else
        J = Jacobi(S)
        LeftIntegralWrapper(LeftIntegral(J,k)*Conversion(S,J),0.5)
    end
end

LeftIntegral(S::Chebyshev,k) = LeftIntegralWrapper(
    LeftIntegral(Ultraspherical(1//2,domain(S)),k)*Conversion(S,Ultraspherical(1//2,domain(S))),
    0.5)


function rangespace{T,DD<:Interval}(Q::ConcreteLeftIntegral{Jacobi{T,DD},Float64})
    μ=Q.order
    S=domainspace(Q)

    JacobiWeight(S.b+μ,0.,Jacobi(S.a-μ,S.b+μ,domain(S)))
end

function RightIntegral(S::Jacobi,k)
    if S.a==0
        ConcreteRightIntegral(S,k)
    else
        J=Jacobi(0.,S.b,domain(S))
        RightIntegralWrapper(RightIntegral(J,k)*Conversion(S,J),k)
    end
end

function rangespace{T,DD<:Interval}(Q::ConcreteRightIntegral{Jacobi{T,DD},Float64})
    μ=Q.order
    S=domainspace(Q)
    @assert S.a==0

    JacobiWeight(0.,S.a+μ,Jacobi(S.a+μ,S.b-μ,domain(S)))
end

for TYP in (:ConcreteLeftIntegral,:ConcreteRightIntegral)
    @eval begin
        bandinds{T,DD<:Interval}(Q::$TYP{Jacobi{T,DD},Float64})=(0,0)
        getindex{T,DD<:Interval}(Q::$TYP{Jacobi{T,DD},Float64},k::Integer,j::Integer) =
            jacobi_frac_getindex(domain(Q),0.,Q.order,k,j)
    end
end


jacobi_frac_getindex(d::Interval,α,μ,k::Integer,j::Integer) =
    jacobi_frac_getindex((arclength(d)/2)^μ,α,μ,k,j)
jacobi_frac_getindex(c::Number,α,μ,k::Integer,j::Integer) =
    k==j ? c*exp(lgamma(α+k)-lgamma(α+μ+k)) : zero(promote_type(typeof(c),typeof(α),typeof(μ)))


# jacobi_frac_addentries!(d::Interval,α,μ,A,kr::UnitRange)=
#     jacobi_frac_addentries!((arclength(d)/2)^μ,α,μ,A,kr)
# function jacobi_frac_addentries!(c::Number,α,μ,A,kr::UnitRange)
#     γ=c*gamma(α+1)/gamma(α+1+μ)
#     for k=1:first(kr)-1
#         γ*=(α+k)/(α+μ+k)
#     end
#     for k=kr
#        A[k,k]+=γ
#        γ*=(α+k)/(α+μ+k)
#         #should be gamma(S.α+k)/gamma(S.α+μ+k)=
#         #(S.α+k-1)/(S.α+μ+k-1)*gamma(S.α+k-1)/gamma(S.α+μ+k-1)
#     end
#     A
# end


function LeftIntegral{T,DD}(S::JacobiWeight{Jacobi{T,DD}},k)
    J=S.space
    @assert S.β==0
    @assert S.α==J.b
    ConcreteLeftIntegral(S,k)
end

function RightIntegral{T,DD}(S::JacobiWeight{Jacobi{T,DD}},k)
    J=S.space
    @assert S.β==J.a
    @assert S.α==0
    ConcreteRightIntegral(S,k)
end



function LeftIntegral{DD}(S::JacobiWeight{Chebyshev{DD}},k)
    # convert to Jacobi
    Q=LeftIntegral(JacobiWeight(S.α,S.β,Jacobi(S.space)),k)
    LeftIntegralWrapper(Q*Conversion(S,domainspace(Q)),k)
end

function RightIntegral{DD}(S::JacobiWeight{Chebyshev{DD}},k)
    # convert to Jacobi
    Q=RightIntegral(JacobiWeight(S.α,S.β,Jacobi(S.space)),k)
    RightIntegralWrapper(Q*Conversion(S,domainspace(Q)),k)
end

for (TYP,WRAP) in ((:LeftIntegral,:LeftIntegralWrapper),
                    (:RightIntegral,:RightIntegralWrapper))
    @eval function $TYP{PS<:PolynomialSpace}(S::JacobiWeight{PS},k)
        JS = JacobiWeight(S.α,S.β,Jacobi(S.space))
        $WRAP($TYP(JS,k)*Conversion(S,JS),k)
    end
end


#DLMF18.17.9
function rangespace{T,DD<:Interval}(Q::ConcreteLeftIntegral{JacobiWeight{Jacobi{T,DD},DD},Float64})
    μ=Q.order
    S=domainspace(Q)
    J=S.space
    if isapprox(S.α,-μ)
        Jacobi(J.a-μ,J.b+μ,domain(J))
    else
        JacobiWeight(S.α+μ,0.,Jacobi(J.a-μ,J.b+μ,domain(J)))
    end
end

function rangespace{T,DD<:Interval}(Q::ConcreteRightIntegral{JacobiWeight{Jacobi{T,DD},DD},Float64})
    μ=Q.order
    S=domainspace(Q)
    J=S.space
    @assert S.β==J.a
    @assert S.α==0
    if isapprox(S.β,-μ)
        Jacobi(J.a+μ,J.b-μ,domain(J))
    else
        JacobiWeight(0.,S.β+μ,Jacobi(J.a+μ,J.b-μ,domain(J)))
    end
end

for TYP in (:ConcreteLeftIntegral,:ConcreteRightIntegral)
    @eval bandinds{T,DD<:Interval}(Q::$TYP{JacobiWeight{Jacobi{T,DD},DD},Float64})=(0,0)
end

getindex{T,DD<:Interval}(Q::ConcreteLeftIntegral{JacobiWeight{Jacobi{T,DD},DD},Float64},k::Integer,j::Integer) =
    jacobi_frac_getindex(domain(Q),domainspace(Q).α,Q.order,k,j)

getindex{T,DD<:Interval}(Q::ConcreteRightIntegral{JacobiWeight{Jacobi{T,DD},DD},Float64},k::Integer,j::Integer) =
    jacobi_frac_getindex(domain(Q),domainspace(Q).β,Q.order,k,j)

function choosedomainspace{T<:Float64}(Q::LeftIntegral{UnsetSpace,T},sp::JacobiWeight)
    @assert sp.α>0 && isapproxinteger(sp.β)

    if isapprox(Q.order,sp.α)
        Legendre(domain(sp))
    else
        JacobiWeight(sp.α-Q.order,0.,Jacobi(0.,sp.α-Q.order,domain(sp)))
    end
end

choosedomainspace{T<:Float64}(Q::LeftIntegral{UnsetSpace,T},sp::PolynomialSpace)=
                JacobiWeight(-Q.order,0.,Jacobi(-Q.order,-Q.order,domain(sp)))



function choosedomainspace{T<:Float64}(Q::RightIntegral{UnsetSpace,T},sp::JacobiWeight)
    @assert sp.β>0  && isapproxinteger(sp.α)

    if isapprox(Q.order,sp.β)
        Legendre(domain(sp))
    else
        JacobiWeight(0.,sp.β-Q.order,Jacobi(sp.β-Q.order,0.))
    end
end
choosedomainspace{T<:Float64}(Q::RightIntegral{UnsetSpace,T},sp::PolynomialSpace)=
    JacobiWeight(0.,-Q.order,Jacobi(-Q.order,-Q.order,domain(sp)))




# Define Left/RightDerivative
for (DTYP,QTYP,DWRAP,QWRAP) in ((:LeftDerivative,:LeftIntegral,:LeftDerivativeWrapper,:LeftIntegralWrapper),
                            (:RightDerivative,:RightIntegral,:RightDerivativeWrapper,:RightIntegralWrapper))
    @eval begin
        function $DTYP(S::Space,k::Real)
            i=ceil(Int,k)
            r=i-k
            $DWRAP(i<0?$QTYP(S,-k):Derivative(i)*$QTYP(S,r),k)
        end
        $QTYP(S::SumSpace,k) = $QWRAP(InterlaceOperator(Diagonal([map(s->$QTYP(s,k),S.spaces)...]),SumSpace),k)
    end
end
