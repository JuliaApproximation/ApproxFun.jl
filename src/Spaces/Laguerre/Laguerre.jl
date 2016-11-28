export Laguerre, LaguerreWeight, WeightedLaguerre

#####
# recα/β/γ are given by
#       x p_{n-1} =γ_n p_{n-2} + α_n p_{n-1} +  p_n β_n
# x L_n^{(α)}(x)  = - (n+α) L_{n-1}^{(α)}(x) +
#                    (2 n + α+ 1) L_n^{(α)}(x)-(n+1)L_{n+1}^{(α)}(x)
# x L_{n-1}^{(α)}(x)  = - (n-1+α) L_{n-2}^{(α)}(x) +
#                    (2 n + α- 1) L_{n-1}^{(α)}(x)-n L_{n}^{(α)}(x)
#####

#####
# jacobirecA/B/C is from dlmf:
# p_{n+1} = (A_n x + B_n)p_n - C_n p_{n-1}
#####

immutable Laguerre{T} <: PolynomialSpace{Ray{false,Float64}}
    α::T
end

Laguerre() = Laguerre(0)

spacescompatible(A::Laguerre,B::Laguerre) = A.α ≈ B.α

canonicaldomain(::Laguerre) = Ray()
domain(::Laguerre) = Ray()
tocanonical(::Laguerre,x) = x

recα{T}(::Type{T},L::Laguerre,k) = T(2k+L.α-1)
recβ{T}(::Type{T},L::Laguerre,k) = T(-k)
recγ{T}(::Type{T},L::Laguerre,k) = T(-(k-1+L.α))


recA{T}(::Type{T},L::Laguerre,k) = T(-1/(k+1))
recB{T}(::Type{T},L::Laguerre,k) = T((2k+L.α+1)/(k+1))
recC{T}(::Type{T},L::Laguerre,k) = T((k+L.α)/(k+1))



Derivative(L::Laguerre,k) =
    DerivativeWrapper(SpaceOperator(ToeplitzOperator(Float64[],[zeros(k);(-1.)^k]),
                                    L,Laguerre(L.α+k)))


union_rule(A::Laguerre,B::Laguerre) = Laguerre(min(A.α,B.α))
maxspace_rule(A::Jacobi,B::Jacobi) = Laguerre(max(A.α,B.α))
function conversion_rule(A::Laguerre,B::Laguerre)
    if !isapproxinteger(A.α-B.α)
        NoSpace()
    else
        Laguerre(min(A.α,B.α))
    end
end




function Conversion(A::Laguerre,B::Laguerre)
    @assert isapproxinteger(A.α - B.α)
    @assert B.α > A.α
    if B.α == A.α+1
        ConversionWrapper(SpaceOperator(
                ToeplitzOperator(Float64[],[1.,-1.]),
                            A,B))
    else
        Conversion(A,Laguerre(A.α+1),B)
    end
end


# x^α*exp(-L*x)
immutable LaguerreWeight{S,T} <: WeightSpace{S,RealBasis,Ray{false,Float64},1}
    space::S
    α::T
    L::T
end

LaguerreWeight(space,α) = LaguerreWeight(space,α,1.0)

WeightedLaguerre(α) = LaguerreWeight(Laguerre(α),α)

spacescompatible(a::LaguerreWeight,b::LaguerreWeight) =
    spacescompatible(a.space,b.space) && isapprox(a.α,b.α) && isapprox(a.L,b.L)

function Derivative(sp::LaguerreWeight,k)
    @assert sp.α == 0
    if k==1
        D=Derivative(sp.space)
        D2=D-sp.L*I
        DerivativeWrapper(SpaceOperator(D2,sp,LaguerreWeight(rangespace(D2),sp.α,sp.L)),1)
    else
        D=Derivative(sp)
        DerivativeWrapper(TimesOperator(Derivative(rangespace(D),k-1).op,D.op),k)
    end
end

weight(L::LaguerreWeight,x) = x.^L.α.*exp(-L.L*x)
