export Hermite,GaussWeight

#TODO: Add general lines
immutable Hermite{T} <: PolynomialSpace{Line{T}}
    L::T
end
Hermite()=Hermite(1.0)

domain(::Hermite)=Line()
canonicalspace(H::Hermite)=Hermite()
spacescompatible(::Hermite,::Hermite)=true #TODO:L



#####
# recα/β/γ are given by
#       x p_{n-1} =γ_n p_{n-2} + α_n p_{n-1} +  p_n β_n
#####


recα(::Type,::Hermite,k)=0;recβ(::Type,::Hermite,k)=0.5;recγ(::Type,::Hermite,k)=k-1

bandinds{H<:Hermite}(D::Derivative{H})=0,D.order
rangespace{H<:Hermite}(D::Derivative{H})=domainspace(D)
function addentries!{H<:Hermite}(D::Derivative{H},A,kr::Range)
    m = D.order
    C = 2^m
    for k=kr
        A[k,k+D.order]+=C*pochhammer(k,m)
    end
    A
end

function hermitep(r::Range,x::Number)
    n=r[end]+1
    T = typeof(x)
    H = Hermite()
    if n≤2
        v=[1.,2x]
    else
        v=Array(promote_type(Float64,typeof(x)),n)  # x may be complex
        v[1]=1.
        v[2]=2x

        for k=2:n-1
            v[k+1]=((x-recα(T,H,k))*v[k] - recγ(T,H,k)*v[k-1])/recβ(T,H,k)
        end
    end
    v[r+1]
end
hermitep(n::Integer,v::Number)=hermitep(n:n,v)[1]
hermitep(n::Range,v::Vector)=transpose(hcat(map(x->hermitep(n,x),v)...))
hermitep(n::Integer,v::Vector)=map(x->hermitep(n,x),v)


identity_fun(sp::Hermite)=Fun([0.,0.5],sp)


# exp(-Lx^2)
immutable GaussWeight{S,T} <: WeightSpace
    space::S
    L::T
end
spacescompatible(a::GaussWeight,b::GaussWeight)=spacescompatible(a.space,b.space)&&isapprox(a.L,b.L)

function Derivative(sp::GaussWeight,k)
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

weight(H::GaussWeight,x)=exp(-H.L*x.^2)

include("hermitetransform.jl")
