export Jacobi,Legendre,WeightedJacobi


immutable Jacobi{T,D<:Domain} <: PolynomialSpace{D}
    a::T
    b::T
    domain::D
end
Legendre(domain)=Jacobi(0.,0.,domain)
Legendre()=Legendre(Interval())
Jacobi(a,b,d::Domain)=Jacobi{promote_type(typeof(a),typeof(b)),typeof(d)}(a,b,d)
Jacobi(a,b,d)=Jacobi(a,b,Domain(d))
Jacobi(a,b)=Jacobi(a,b,Interval())
Jacobi{m}(A::Ultraspherical{m})=Jacobi(m-0.5,m-0.5,domain(A))


Base.promote_rule{T,V,D}(::Type{Jacobi{T,D}},::Type{Jacobi{V,D}})=Jacobi{promote_type(T,V),D}
Base.convert{T,V,D}(::Type{Jacobi{T,D}},J::Jacobi{V,D})=Jacobi{T,D}(J.a,J.b,J.domain)

typealias WeightedJacobi{D} JacobiWeight{Jacobi{Float64,D},D}

Base.call(::Type{WeightedJacobi},α,β,d::Domain)=JacobiWeight(α,β,Jacobi(β,α,d))
Base.call(::Type{WeightedJacobi},α,β)=JacobiWeight(α,β,Jacobi(β,α))

spacescompatible(a::Jacobi,b::Jacobi)=a.a==b.a && a.b==b.b

function canonicalspace(S::Jacobi)
    if isapproxinteger(S.a+0.5) && isapproxinteger(S.b+0.5)
        Chebyshev(domain(S))
    else
        # return space with parameters in (-1,0.]
        Jacobi(mod(S.a,-1),mod(S.b,-1),domain(S))
    end
end

#####
# jacobirecA/B/C is from dlmf:
# p_{n+1} = (A_n x + B_n)p_n - C_n p_{n-1}
#####
jacobirecA{T}(::Type{T},α,β,k)=k==0&&((α+β==0)||(α+β==-1))?.5*(α+β)+one(T):(2k+α+β+one(T))*(2k+α+β+2one(T))/(2*(k+one(T))*(k+α+β+one(T)))
jacobirecB{T}(::Type{T},α,β,k)=k==0&&((α+β==0)||(α+β==-1))?.5*(α-β)*one(T):(α-β)*(α+β)*(2k+α+β+one(T))/(2*(k+one(T))*(k+α+β+one(T))*(2one(T)*k+α+β))
jacobirecC{T}(::Type{T},α,β,k)=(one(T)*k+α)*(one(T)*k+β)*(2k+α+β+2one(T))/((k+one(T))*(k+α+β+one(T))*(2one(T)*k+α+β))

#####
# jacobirecA/B/C is from dlmf:
# x p_{n-1} =γ_n p_{n-2} + α_n p_{n-1} +  p_n β_n
#####

jacobirecγ{T}(::Type{T},α,β,k)=jacobirecC(T,α,β,k-1)/jacobirecA(T,α,β,k-1)
jacobirecα{T}(::Type{T},α,β,k)=-jacobirecB(T,α,β,k-1)/jacobirecA(T,α,β,k-1)
jacobirecβ{T}(::Type{T},α,β,k)=1/jacobirecA(T,α,β,k-1)

for (REC,JREC) in ((:recα,:jacobirecα),(:recβ,:jacobirecβ),(:recγ,:jacobirecγ),
                   (:recA,:jacobirecA),(:recB,:jacobirecB),(:recC,:jacobirecC))
    @eval $REC{T}(::Type{T},sp::Jacobi,k)=$JREC(T,sp.a,sp.b,k)  #TODO: implement typing
end


function jacobip(r::Range,α,β,x)
    n=r[end]+1
    if n<=2
        v=[1.,.5*(α-β+(2+α+β)*x)]
    else
        T=promote_type(Float64,typeof(x))
        v=Vector{T}(n)  # x may be complex
        v[1]=1.
        v[2]=.5*(α-β+(2+α+β)*x)

        for k=2:n-1
            v[k+1]=((x-jacobirecα(T,α,β,k))*v[k] - jacobirecγ(T,α,β,k)*v[k-1])/jacobirecβ(T,α,β,k)
        end
    end
    v[r+1]
end


function jacobip(r::Range,α,β,x::Number)
    if x==1. && α==0.
        ones(length(r))
    elseif x==-1. && β==0.
        (-1.).^r
    else
        n=r[end]+1
        if n<=2
            v=[1.,.5*(α-β+(2+α+β)*x)]
        else
            T=promote_type(Float64,typeof(x))
            v=Vector{T}(n)  # x may be complex
            v[1]=1.
            v[2]=.5*(α-β+(2+α+β)*x)

            for k=2:n-1
                v[k+1]=((x-jacobirecα(T,α,β,k))*v[k] - jacobirecγ(T,α,β,k)*v[k-1])/jacobirecβ(T,α,β,k)
            end
        end
        v[r+1]
    end
end
jacobip(n::Integer,α,β,v)=jacobip(n:n,α,β,v)[1]
jacobip(n::Range,α,β,v::Vector)=transpose(hcat(map(x->jacobip(n,α,β,x),v)...))
jacobip(n::Integer,α,β,v::Vector)=map(x->jacobip(n,α,β,x),v)
jacobip(n,S::Jacobi,v)=jacobip(n,S.a,S.b,v)





include("jacobitransform.jl")
include("JacobiOperators.jl")



for op in (:(Base.ones),:(Base.zeros))
    @eval ($op){T<:Number}(::Type{T},S::Jacobi)=Fun(($op)(T,1),S)
    @eval ($op)(S::Jacobi)=Fun(($op)(1),S)
end

function identity_fun(J::Jacobi)
    if domain(J)==Interval()
        Fun([(J.b-J.a)/(2+J.a+J.b),2.0/(2+J.a+J.b)],J)
    else
        d=domain(J)
        complexlength(d)/2*(Fun([(J.b-J.a)/(2+J.a+J.b),2.0/(2+J.a+J.b)],J)+1.)+first(d)
    end
end


setdomain(S::Jacobi,d::Domain)=Jacobi(S.a,S.b,d)


# O(min(m,n)) Jacobi conjugated inner product

function conjugatedinnerproduct{S,V}(sp::Jacobi,u::Vector{S},v::Vector{V})
    T,mn = promote_type(S,V),min(length(u),length(v))
    α,β = sp.a,sp.b
    if mn > 1
        wi = 2^(α+β+1)*gamma(α+1)*gamma(β+1)/gamma(α+β+2)
        ret = u[1]*wi*v[1]
        for i=2:mn
            wi *= (α+i-1)*(β+i-1)/(i-1)/(i-1+α+β)*(2i+α+β-3)/(2i+α+β-1)
            ret += u[i]*wi*v[i]
        end
        return ret
    elseif mn > 0
        wi = 2^(α+β+1)*gamma(α+1)*gamma(β+1)/gamma(α+β+2)
        return u[1]*wi*v[1]
    else
        return zero(promote_type(eltype(u),eltype(v)))
    end
end

function bilinearform{J<:Jacobi}(f::Fun{J},g::Fun{J})
    @assert domain(f) == domain(g)
    if f.space.a == g.space.a == 0. && f.space.b == g.space.b == 0.
        return complexlength(domain(f))/2*conjugatedinnerproduct(g.space,f.coefficients,g.coefficients)
    else
        return defaultbilinearform(f,g)
    end
end

function bilinearform{J<:Jacobi,DD<:Interval}(f::Fun{JacobiWeight{J,DD}},g::Fun{J})
    @assert domain(f) == domain(g)
    if f.space.β == f.space.space.a == g.space.a && f.space.α == f.space.space.b == g.space.b
        return complexlength(domain(f))/2*conjugatedinnerproduct(g.space,f.coefficients,g.coefficients)
    else
        return defaultbilinearform(f,g)
    end
end

function bilinearform{J<:Jacobi,DD<:Interval}(f::Fun{J},g::Fun{JacobiWeight{J,DD}})
    @assert domain(f) == domain(g)
    if g.space.β == g.space.space.a == f.space.a && g.space.α == g.space.space.b == f.space.b
        return complexlength(domain(f))/2*conjugatedinnerproduct(f.space,f.coefficients,g.coefficients)
    else
        return defaultbilinearform(f,g)
    end
end

function bilinearform{J<:Jacobi,DD<:Interval}(f::Fun{JacobiWeight{J,DD}},g::Fun{JacobiWeight{J,DD}})
    @assert domain(f) == domain(g)
    if f.space.β + g.space.β == f.space.space.a == g.space.space.a && f.space.α + g.space.α == f.space.space.b == g.space.space.b
        return complexlength(domain(f))/2*conjugatedinnerproduct(f.space.space,f.coefficients,g.coefficients)
    else
        return defaultbilinearform(f,g)
    end
end

function linebilinearform{J<:Jacobi}(f::Fun{J},g::Fun{J})
    @assert domain(f) == domain(g)
    if f.space.a == g.space.a == 0. && f.space.b == g.space.b == 0.
        return length(domain(f))/2*conjugatedinnerproduct(g.space,f.coefficients,g.coefficients)
    else
        return defaultlinebilinearform(f,g)
    end
end

function linebilinearform{J<:Jacobi,DD<:Interval}(f::Fun{JacobiWeight{J,DD}},g::Fun{J})
    @assert domain(f) == domain(g)
    if f.space.β == f.space.space.a == g.space.a && f.space.α == f.space.space.b == g.space.b
        return length(domain(f))/2*conjugatedinnerproduct(g.space,f.coefficients,g.coefficients)
    else
        return defaultlinebilinearform(f,g)
    end
end

function linebilinearform{J<:Jacobi,DD<:Interval}(f::Fun{J},g::Fun{JacobiWeight{J,DD}})
    @assert domain(f) == domain(g)
    if g.space.β == g.space.space.a == f.space.a && g.space.α == g.space.space.b == f.space.b
        return length(domain(f))/2*conjugatedinnerproduct(f.space,f.coefficients,g.coefficients)
    else
        return defaultlinebilinearform(f,g)
    end
end

function linebilinearform{J<:Jacobi,DD<:Interval}(f::Fun{JacobiWeight{J,DD}},g::Fun{JacobiWeight{J,DD}})
    @assert domain(f) == domain(g)
    if f.space.β + g.space.β == f.space.space.a == g.space.space.a && f.space.α + g.space.α == f.space.space.b == g.space.space.b
        return length(domain(f))/2*conjugatedinnerproduct(f.space.space,f.coefficients,g.coefficients)
    else
        return defaultlinebilinearform(f,g)
    end
end
