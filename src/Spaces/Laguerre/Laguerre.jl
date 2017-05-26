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
# laguerrerecA/B/C is from dlmf:
# p_{n+1} = (A_n x + B_n)p_n - C_n p_{n-1}
#####

struct Laguerre{T} <: PolynomialSpace{Ray{false,Float64},Float64}
    α::T
end

Laguerre() = Laguerre(0)

spacescompatible(A::Laguerre,B::Laguerre) = A.α ≈ B.α

canonicaldomain(::Laguerre) = Ray()
domain(::Laguerre) = Ray()
tocanonical(::Laguerre,x) = x

@inline laguerrerecα{T}(::Type{T},α,k) = T(2k+α-1)
@inline laguerrerecβ{T}(::Type{T},::,k) = T(-k)
@inline laguerrerecγ{T}(::Type{T},α,k) = T(-(k-1+α))


@inline laguerrerecA{T}(::Type{T},::,k) = T(-1/(k+1))
@inline laguerrerecB{T}(::Type{T},α,k) = T((2k+α+1)/(k+1))
@inline laguerrerecC{T}(::Type{T},α,k) = T((k+α)/(k+1))

for (REC,JREC) in ((:recα,:laguerrerecα),(:recβ,:laguerrerecβ),(:recγ,:laguerrerecγ),
                   (:recA,:laguerrerecA),(:recB,:laguerrerecB),(:recC,:laguerrerecC))
    @eval @inline $REC{T}(::Type{T},sp::Laguerre,k) = $JREC(T,sp.α,k)
end


identity_fun(L::Laguerre) = Fun(L,[1+L.α,-1.0])



function laguerrel{T}(::Type{T},r::Range,α,x)
    if isempty(r)
        T[]
    else
        n=r[end]+1
        if n<=2
            v=T[1,1.0-x+α]
        else
            v=Vector{T}(n)  # x may be complex
            v[1]=1
            v[2]=1.0-x+α

            @inbounds for k=2:n-1
                v[k+1]=((x-laguerrerecα(T,α,k))*v[k] - laguerrerecγ(T,α,k)*v[k-1])/laguerrerecβ(T,α,k)
            end
        end
        v[r+1]
    end
end

laguerrel(r::Range,α,x) = laguerrel(promote_type(typeof(α),typeof(x)),r,α,x)

laguerrel{T}(::Type{T},n::Integer,α,v) = laguerrel(T,n:n,α,v)[1]
laguerrel(n::Integer,α,v) = laguerrel(n:n,α,v)[1]
laguerrel{T}(::Type{T},n::Range,α,v::AbstractVector) = transpose(hcat(map(x->laguerrel(T,n,α,x),v)...))
laguerrel(n::Range,α,v::AbstractVector) = transpose(hcat(map(x->laguerrel(n,α,x),v)...))
laguerrel{T}(::Type{T},n::Integer,α,v::AbstractVector) = map(x->laguerrel(T,n,α,x),v)
laguerrel(n::Integer,α,v::AbstractVector) = map(x->laguerrel(n,α,x),v)
laguerrel{T}(::Type{T},n::Integer,S::Laguerre,v::AbstractVector) = laguerrel(T,n,S.a,S.b,v)
laguerrel{T}(::Type{T},n::Range,S::Laguerre,v::AbstractVector) = laguerrel(T,n,S.a,S.b,v)
laguerrel{T}(::Type{T},n,S::Laguerre,v::AbstractVector) = laguerrel(T,n,S.a,S.b,v)
laguerrel{T}(::Type{T},n::Integer,S::Laguerre,v) = laguerrel(T,n,S.a,S.b,v)
laguerrel{T}(::Type{T},n::Range,S::Laguerre,v) = laguerrel(T,n,S.a,S.b,v)
laguerrel{T}(::Type{T},n,S::Laguerre,v) = laguerrel(T,n,S.a,S.b,v)
laguerrel(n::Integer,S::Laguerre,v::AbstractVector) = laguerrel(n,S.a,S.b,v)
laguerrel(n::Range,S::Laguerre,v::AbstractVector) = laguerrel(n,S.a,S.b,v)
laguerrel(n,S::Laguerre,v::AbstractVector) = laguerrel(n,S.a,S.b,v)
laguerrel(n::Integer,S::Laguerre,v) = laguerrel(n,S.a,S.b,v)
laguerrel(n::Range,S::Laguerre,v) = laguerrel(n,S.a,S.b,v)
laguerrel(n,S::Laguerre,v) = laguerrel(n,S.a,S.b,v)


struct LaguerreTransformPlan{T,TT}
    space::Laguerre{TT}
    points::Vector{T}
    weights::Vector{T}
end

plan_transform(S::Laguerre,v::AbstractVector) = LaguerreTransformPlan(S,gausslaguerre(length(v),1.0S.α)...)
function *(plan::LaguerreTransformPlan,vals)
#    @assert S==plan.space
    x,w = plan.points, plan.weights
    V=laguerrel(0:length(vals)-1,plan.space.α,x)'
    #w2=w.*x.^(S.α-plan.space.α)   # need to weight if plan is different
    w2=w
    nrm=(V.^2)*w2
    V*(w2.*vals)./nrm
end


points(L::Laguerre,n) = gausslaguerre(n,1.0L.α)[1]



Derivative(L::Laguerre,k) =
    DerivativeWrapper(SpaceOperator(ToeplitzOperator(Float64[],[zeros(k);(-1.)^k]),
                                    L,Laguerre(L.α+k)))


union_rule(A::Laguerre,B::Laguerre) = Laguerre(min(A.α,B.α))
maxspace_rule(A::Laguerre,B::Laguerre) = Laguerre(max(A.α,B.α))
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
struct LaguerreWeight{S,T} <: WeightSpace{S,Ray{false,Float64},Float64}
    α::T
    L::T
    space::S
end

LaguerreWeight(α,space::Space) = LaguerreWeight(α,1.0,space)

WeightedLaguerre(α) = LaguerreWeight(α,Laguerre(α))

@inline laguerreweight(α,L,x) = x.^α.*exp(-L*x)
@inline weight(L::LaguerreWeight,x) = laguerreweight(L.α,L.L,x)


spacescompatible(a::LaguerreWeight,b::LaguerreWeight) =
    spacescompatible(a.space,b.space) && isapprox(a.α,b.α) && isapprox(a.L,b.L)

function Derivative(sp::LaguerreWeight,k)
    @assert sp.α == 0
    if k==1
        D=Derivative(sp.space)
        D2=D-sp.L*I
        DerivativeWrapper(SpaceOperator(D2,sp,LaguerreWeight(sp.α,sp.L,rangespace(D2))),1)
    else
        D=Derivative(sp)
        DerivativeWrapper(TimesOperator(Derivative(rangespace(D),k-1).op,D.op),k)
    end
end

function Multiplication(f::Fun,S::LaguerreWeight)
    M=Multiplication(f,S.space)
    rsp=LaguerreWeight(S.α,S.L,rangespace(M))
    MultiplicationWrapper(f,SpaceOperator(M,S,rsp))
end


function conversion_rule(A::LaguerreWeight,B::LaguerreWeight)
    if isapproxinteger(A.α-B.α) && A.L == B.L
        ct=conversion_type(A.space,B.space)
        ct==NoSpace()?NoSpace():LaguerreWeight(max(A.α,B.α),A.L,ct)
    else
        NoSpace()
    end
end


conversion_rule{D<:Ray}(A::LaguerreWeight,B::Space{D}) = conversion_type(A,LaguerreWeight(0,0,B))


function Conversion(A::LaguerreWeight,B::LaguerreWeight)
    @assert isapproxinteger(A.α-B.α) && A.L == B.L

    if isapprox(A.α,B.α)
        ConversionWrapper(SpaceOperator(Conversion(A.space,B.space),A,B))
    else
        @assert A.α≥B.α
        # first check if a multiplication by LaguerreWeight times B.space is overloaded
        # this is currently designed for Laguerre multiplied by (1-x), etc.
        αdif=round(Int,A.α-B.α)

        M=Multiplication(laguerreweight(αdif,A.L,d),
                         A.space)

        if rangespace(M) == LaguerreWeight(αdif,A.L,A.space)
            # M is the default, so we should use multiplication by polynomials instead
            x=Fun(identity,A.space)
            m=x^αdif
            MC=promoterangespace(Multiplication(m,A.space),B.space)

            ConversionWrapper(SpaceOperator(MC,A,B))# Wrap the operator with the correct spaces
        else
            ConversionWrapper(SpaceOperator(promoterangespace(M,B.space),
                                            A,B))
        end
    end
end


Conversion{D<:Ray}(A::ConstantSpace{D},B::LaguerreWeight) = error("Cannot convert constants to LaguerreWeight.")
Conversion(a::SubSpace{S,IT,DD,RR},b::S) where {S<:LaguerreWeight,IT,DD<:Ray,RR} =
    ConcreteConversion(a,b)
Conversion{D<:Ray,RR<:Real}(A::Space{D,RR},B::LaguerreWeight) = ConversionWrapper(
    SpaceOperator(
        Conversion(LaguerreWeight(0,0,A),B),
        A,B))
Conversion{D<:Ray,RR<:Real}(A::LaguerreWeight,B::Space{D,RR}) = ConversionWrapper(
    SpaceOperator(
        Conversion(A,LaguerreWeight(0,0,B)),
        A,B))
