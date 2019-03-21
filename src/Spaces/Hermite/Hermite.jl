export Hermite,NormalizedHermite,GaussWeight

#TODO: Add general lines

"""
`Hemite(L)` represents `H_k(sqrt(L) * x)` where `H_k` are Hermite polynomials.
`Hermite()` is equivalent to `Hermite(1.0)`.
"""
struct Hermite{T} <: PolynomialSpace{Line{false,Float64},Float64}
    L::T
end
Hermite() = Hermite(1.0)

NormalizedHermite() = NormalizedPolynomialSpace(Hermite())
NormalizedHermite(L) = NormalizedPolynomialSpace(Hermite(L))

domain(::Hermite) = Line()
canonicalspace(H::Hermite) = Hermite()
spacescompatible(::Hermite,::Hermite) = true #TODO:L
canonicaldomain(H::Hermite) = Line()
tocanonical(H::Hermite,x) = x


#####
# recα/β/γ are given by
#       x p_{n-1} =γ_n p_{n-2} + α_n p_{n-1} +  p_n β_n
#####


recα(::Type,::Hermite,k) = 0;
recβ(::Type,H::Hermite,k) = 0.5 / sqrt(H.L);
recγ(::Type,H::Hermite,k) = (k-1) / sqrt(H.L);
recA(::Type,H::Hermite,k) = 2 * sqrt(H.L);
recB(::Type,::Hermite,k) = 0;
recC(::Type,::Hermite,k) = 2k;

normalization(::Type{T}, sp::Hermite, k::Int) where T = (warn("This normalization would have overflown sooner than you wished! Not normalizing."); one(T))

Derivative(H::Hermite,order) = ConcreteDerivative(H,order)


bandwidths(D::ConcreteDerivative{H}) where {H<:Hermite} = 0,D.order
rangespace(D::ConcreteDerivative{H}) where {H<:Hermite} = domainspace(D)
getindex(D::ConcreteDerivative{H},k::Integer,j::Integer) where {H<:Hermite} =
        j==k+D.order ? one(eltype(D))*2^D.order*pochhammer(k,D.order) : zero(eltype(D))

function hermitep(r::AbstractRange,x::Number)
    n = r[end] + 1

    T = typeof(x)
    H = Hermite()
    if n ≤ 2
        v = [1.,2x]
    else
        v=Array{promote_type(Float64,typeof(x))}(undef, n)  # x may be complex
        v[1] = 1.
        v[2] = 2x

        for k = 2:n-1
            v[k + 1] = ((x - recα(T,H,k)) * v[k] - recγ(T,H,k) * v[k - 1]) / recβ(T,H,k)
        end
    end
    v[r.+1]
end
hermitep(n::Integer, v::Number) = hermitep(n:n,v)[1]




# exp(-Lx^2)
struct GaussWeight{S,T} <: WeightSpace{S,Line{Float64},Float64}
    space::S
    L::T
end

GaussWeight(H::Hermite) = GaussWeight(H,H.L)
GaussWeight() = GaussWeight(Hermite())

"""
`GaussWeight(Hermite(L), L)` is a space spanned by `exp(-Lx²) * H_k(sqrt(L) * x)`
where `H_k(x)`'s are Hermite polynomials.

`GaussWeight()` is equivalent to `GaussWeight(Hermite(), 1.0)` by default.
"""
Fun(::typeof(identity), sp::Hermite) = Fun(sp,[0.,0.5])
Fun(::typeof(identity), sp::GaussWeight) = Fun(identity, sp.space)

spacescompatible(a::GaussWeight,b::GaussWeight)=spacescompatible(a.space,b.space)&&isapprox(a.L,b.L)

function Derivative(sp::GaussWeight,k::Integer)
   if k == 1
        x = Multiplication(Fun(identity,sp.space),sp.space)
        D = Derivative(sp.space)
        D2 = D - (2sp.L)x
        DerivativeWrapper(SpaceOperator(D2,sp,GaussWeight(rangespace(D2),sp.L)),1)
    else
        D = Derivative(sp)
        DerivativeWrapper(TimesOperator(Derivative(rangespace(D),k-1).op,D.op),k)
    end
end

weight(H::GaussWeight,x) = exp(-H.L * x^2)

function Base.sum(f::Fun{GaussWeight{H,T}}) where {H<:Hermite,T}
    @assert space(f).space.L==space(f).L  # only implemented with matching weight
    f.coefficients[1]*sqrt(convert(T,π))/sqrt(space(f).L)
end

include("hermitetransform.jl")





function Multiplication(f::Fun{H}, S::GaussWeight{H}) where H<:Hermite
    M = Multiplication(f, S.space)
    rs = rangespace(M)
    MultiplicationWrapper(f, SpaceOperator(M, S, GaussWeight(rs, S.L)))
end

function Multiplication(f::Fun{GaussWeight{H,T}}, S::Hermite) where {H<:Hermite,T}
    M = Multiplication(Fun(space(f).space, f.coefficients), S)
    rs = rangespace(M)
    MultiplicationWrapper(f, SpaceOperator(M, S, GaussWeight(rs, space(f).L)))
end





function integrate(f::Fun{GW}) where GW<:GaussWeight{H} where H<:Hermite
    n = length(f.coefficients);
    L = space(f).L;
    if L == 0
        return Fun(Hermite(space(f).space.L), [0; f.coefficients] ./ [1; 2:2:2n])
    elseif space(f).space.L != L
        throw(ArgumentError("`integrate` is applicable if only parameters of GaussWeight and Hermite are equal."))
    else
        if n == 0
            return Fun(0)
        elseif f.coefficients[1] == 0
            return Fun(GaussWeight(Hermite(L), L), f.coefficients[2:end] / -sqrt(L))
        else
            g = Fun(GaussWeight(Hermite(L), L), f.coefficients[2:end] / -sqrt(L));
            f₀ = Fun(GaussWeight(Hermite(L), L), [f.coefficients[1]]);
            f₀ = Fun(f₀, Chebyshev(Line()));
            g₀ = integrate(f₀);
            g₀ = g₀ - last(g₀);
            return g + g₀
        end
    end
end
