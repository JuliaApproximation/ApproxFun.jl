export Hermite,GaussWeight

#TODO: Add general lines
struct Hermite{T} <: PolynomialSpace{Line{false,Float64},Float64}
    L::T
end
Hermite()=Hermite(1.0)

domain(::Hermite) = Line()
canonicalspace(H::Hermite) = Hermite()
spacescompatible(::Hermite,::Hermite)=true #TODO:L
canonicaldomain(H::Hermite)=Line()
tocanonical(H::Hermite,x)=x


#####
# recα/β/γ are given by
#       x p_{n-1} =γ_n p_{n-2} + α_n p_{n-1} +  p_n β_n
#####


recα(::Type,::Hermite,k)=0;recβ(::Type,::Hermite,k)=0.5;recγ(::Type,::Hermite,k)=k-1
recA(::Type,::Hermite,k)=2;recB(::Type,::Hermite,k)=0;recC(::Type,::Hermite,k)=2k


Derivative(H::Hermite,order)=ConcreteDerivative(H,order)


bandinds(D::ConcreteDerivative{H}) where {H<:Hermite}=0,D.order
rangespace(D::ConcreteDerivative{H}) where {H<:Hermite}=domainspace(D)
getindex(D::ConcreteDerivative{H},k::Integer,j::Integer) where {H<:Hermite} =
        j==k+D.order?one(eltype(D))*2^D.order*pochhammer(k,D.order):zero(eltype(D))



function hermitep(r::Range,x::Number)
    n=r[end]+1
    T = typeof(x)
    H = Hermite()
    if n≤2
        v=[1.,2x]
    else
        v=Array{promote_type(Float64,typeof(x))}(n)  # x may be complex
        v[1]=1.
        v[2]=2x

        for k=2:n-1
            v[k+1]=((x-recα(T,H,k))*v[k] - recγ(T,H,k)*v[k-1])/recβ(T,H,k)
        end
    end
    v[r+1]
end
hermitep(n::Integer,v::Number) = hermitep(n:n,v)[1]

Fun(::typeof(identity), sp::Hermite) = Fun(sp,[0.,0.5])


# exp(-Lx^2)
struct GaussWeight{S,T} <: WeightSpace{S,Line{Float64},Float64}
    space::S
    L::T
end

GaussWeight(H::Hermite)=GaussWeight(H,H.L)
GaussWeight()=GaussWeight(Hermite())


Fun(::typeof(identity), sp::GaussWeight) = Fun(identity, sp.space)

spacescompatible(a::GaussWeight,b::GaussWeight)=spacescompatible(a.space,b.space)&&isapprox(a.L,b.L)

function Derivative(sp::GaussWeight,k::Integer)
   if k==1
        x=Multiplication(Fun(identity,sp.space),sp.space)
        D=Derivative(sp.space)
        D2=D-(2sp.L)x
        DerivativeWrapper(SpaceOperator(D2,sp,GaussWeight(rangespace(D2),sp.L)),1)
    else
        D=Derivative(sp)
        DerivativeWrapper(TimesOperator(Derivative(rangespace(D),k-1).op,D.op),k)
    end
end

weight(H::GaussWeight,x) = exp(-H.L*x^2)

function Base.sum(f::Fun{GaussWeight{H,T}}) where {H<:Hermite,T}
    @assert space(f).space.L==space(f).L  # only implemented with matching weight
    f.coefficients[1]*sqrt(π)/sqrt(space(f).L)
end

include("hermitetransform.jl")





function Multiplication(f::Fun{H},S::GaussWeight{H}) where H<:Hermite
    M=Multiplication(f,S.space)
    rs=rangespace(M)
    MultiplicationWrapper(f,SpaceOperator(M,S,GaussWeight(rs,rs.L)))
end

function Multiplication(f::Fun{GaussWeight{H,T}},S::Hermite) where {H<:Hermite,T}
    M=Multiplication(Fun(space(f).space,f.coefficients),S)
    rs=rangespace(M)
    MultiplicationWrapper(f,SpaceOperator(M,S,GaussWeight(rs,rs.L)))
end
