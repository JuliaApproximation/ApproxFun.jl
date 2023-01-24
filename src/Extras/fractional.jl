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
        J=Jacobi(zero(S.a),S.a,domain(S))
        LeftIntegralWrapper(LeftIntegral(J,k)*Conversion(S,J),k)
    end
end

function LeftIntegral(S::Ultraspherical,k)
    J = Jacobi(S)
    LeftIntegralWrapper(LeftIntegral(J,k)*Conversion(S,J),0.5)
end

LeftIntegral(S::Chebyshev,k) = LeftIntegralWrapper(
    LeftIntegral(Ultraspherical(1//2,domain(S)),k)*Conversion(S,Ultraspherical(1//2,domain(S))),
    0.5)


function rangespace(Q::ConcreteLeftIntegral{<:Jacobi{<:IntervalOrSegment}})
    μ=Q.order
    S=domainspace(Q)

    JacobiWeight(S.b+μ,zero(S.b),Jacobi(S.b+μ,S.a-μ,domain(S)))
end

function RightIntegral(S::Jacobi,k)
    if S.a==0
        ConcreteRightIntegral(S,k)
    else
        J=Jacobi(S.b,zero(S.b),domain(S))
        RightIntegralWrapper(RightIntegral(J,k)*Conversion(S,J),k)
    end
end

function rangespace(Q::ConcreteRightIntegral{<:Jacobi{<:IntervalOrSegment}})
    μ=Q.order
    S=domainspace(Q)
    @assert S.a==0

    JacobiWeight(zero(S.a),S.a+μ,Jacobi(S.b-μ,S.a+μ,domain(S)))
end

for TYP in (:ConcreteLeftIntegral,:ConcreteRightIntegral)
    @eval begin
        bandwidths(Q::$TYP{<:Jacobi{<:IntervalOrSegment}}) = (0,0)
        getindex(Q::$TYP{<:Jacobi{<:IntervalOrSegment},<:Any,T},k::Integer,j::Integer) where {T} =
            convert(T, jacobi_frac_getindex(domain(Q),zero(T),Q.order,k,j))::T
    end
end


jacobi_frac_getindex(d::IntervalOrSegment,α,μ,k::Integer,j::Integer) =
    jacobi_frac_getindex((arclength(d)/2)^μ,α,μ,k,j)

jacobi_frac_getindex(c::Number,α,μ,k::Integer,j::Integer) =
    k==j ? c*exp(loggamma(α+k)-loggamma(α+μ+k)) : zero(promote_type(typeof(c),typeof(α),typeof(μ)))


# jacobi_frac_addentries!(d::IntervalOrSegment,α,μ,A,kr::UnitRange)=
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


function LeftIntegral(S::JacobiWeight{<:Jacobi}, k)
    J=S.space
    @assert S.α ≈ 0
    @assert S.β ≈ J.b
    ConcreteLeftIntegral(S,k)
end

function RightIntegral(S::JacobiWeight{<:Jacobi}, k)
    J=S.space
    @assert S.α ≈ J.a
    @assert S.β ≈ 0
    ConcreteRightIntegral(S,k)
end



function LeftIntegral(S::JacobiWeight{<:Chebyshev},k)
    # convert to Jacobi
    Q=LeftIntegral(JacobiWeight(S.β,S.α,Jacobi(S.space)),k)
    LeftIntegralWrapper(Q*Conversion(S,domainspace(Q)),k)
end

function RightIntegral(S::JacobiWeight{<:Chebyshev},k)
    # convert to Jacobi
    Q=RightIntegral(JacobiWeight(S.β,S.α,Jacobi(S.space)),k)
    RightIntegralWrapper(Q*Conversion(S,domainspace(Q)),k)
end

for (TYP,WRAP) in ((:LeftIntegral,:LeftIntegralWrapper),
                    (:RightIntegral,:RightIntegralWrapper))

    @eval function $TYP(S::JacobiWeight{<:PolynomialSpace}, k)
        JS = JacobiWeight(S.β,S.α,Jacobi(S.space))
        $WRAP($TYP(JS,k)*Conversion(S,JS),k)
    end
end


#DLMF18.17.9
function rangespace(Q::ConcreteLeftIntegral{<:JacobiWeight{<:Jacobi{DD,RR},DD,RR}}) where {DD<:IntervalOrSegment,RR}
    μ=Q.order
    S=domainspace(Q)
    J=S.space
    if isapprox(S.β,-μ)
        Jacobi(J.b+μ,J.a-μ,domain(J))
    else
        JacobiWeight(S.β+μ,zero(S.β),Jacobi(J.b+μ,J.a-μ,domain(J)))
    end
end

function rangespace(Q::ConcreteRightIntegral{<:JacobiWeight{<:Jacobi{DD,RR},DD,RR}}) where {DD<:IntervalOrSegment,RR}
    μ=Q.order
    S=domainspace(Q)
    J=S.space
    @assert S.α==J.a
    @assert S.β==0
    if isapprox(S.β,-μ)
        Jacobi(J.b-μ,J.a+μ,domain(J))
    else
        JacobiWeight(zero(S.α),S.α+μ,Jacobi(J.b-μ,J.a+μ,domain(J)))
    end
end

for TYP in (:ConcreteLeftIntegral,:ConcreteRightIntegral)
    @eval bandwidths(Q::$TYP{<:JacobiWeight{<:Jacobi{DD,RR},DD,RR}}) where {DD<:IntervalOrSegment,RR} = (0,0)
end

getindex(Q::ConcreteLeftIntegral{<:JacobiWeight{<:Jacobi{DD,RR},DD,RR},<:Any,T},k::Integer,j::Integer) where {DD<:IntervalOrSegment,RR,T} =
    convert(T,jacobi_frac_getindex(domain(Q),domainspace(Q).β,Q.order,k,j))::T

getindex(Q::ConcreteRightIntegral{<:JacobiWeight{<:Jacobi{DD,RR},DD,RR},<:Any,T},k::Integer,j::Integer) where {DD<:IntervalOrSegment,RR,T} =
    convert(T,jacobi_frac_getindex(domain(Q),domainspace(Q).α,Q.order,k,j))::T

function choosedomainspace(Q::LeftIntegral{UnsetSpace}, sp::JacobiWeight)
    @assert sp.β>0 && isapproxinteger(sp.α)

    if isapprox(Q.order,sp.β)
        Legendre(domain(sp))
    else
        JacobiWeight(sp.β-Q.order,0.,Jacobi(sp.β-Q.order,0.,domain(sp)))
    end
end

choosedomainspace(Q::LeftIntegral{UnsetSpace,T}, sp::PolynomialSpace) where {T} =
                JacobiWeight(-Q.order,zero(T),Jacobi(-Q.order,-Q.order,domain(sp)))



function choosedomainspace(Q::RightIntegral{UnsetSpace,T},sp::JacobiWeight) where T
    @assert sp.α>0  && isapproxinteger(sp.β)

    if isapprox(Q.order,sp.α)
        Legendre(domain(sp))
    else
        JacobiWeight(zero(T),sp.α-Q.order,Jacobi(0.,sp.α-Q.order))
    end
end
choosedomainspace(Q::RightIntegral{UnsetSpace,T},sp::PolynomialSpace) where {T}=
    JacobiWeight(zero(T),-Q.order,Jacobi(-Q.order,-Q.order,domain(sp)))




# Define Left/RightDerivative
for (DTYP,QTYP,DWRAP,QWRAP) in ((:LeftDerivative,:LeftIntegral,:LeftDerivativeWrapper,:LeftIntegralWrapper),
                            (:RightDerivative,:RightIntegral,:RightDerivativeWrapper,:RightIntegralWrapper))
    @eval begin
        function $DTYP(S::Space,k::Real)
            i=ceil(Int,k)
            r=i-k
            $DWRAP(i<0 ? $QTYP(S,-k) : Derivative(i)*$QTYP(S,r),k)
        end
        $QTYP(S::SumSpace,k) = $QWRAP(InterlaceOperator(Diagonal([map(s->$QTYP(s,k),S.spaces)...]),SumSpace),k)
    end
end
