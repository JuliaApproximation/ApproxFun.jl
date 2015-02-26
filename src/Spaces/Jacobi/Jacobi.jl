export Jacobi,Legendre

immutable Jacobi <: PolynomialSpace
    a::Float64
    b::Float64
    domain::Union(IntervalDomain,AnyDomain)
end
Legendre(domain)=Jacobi(0.,0.,domain)
Legendre()=Legendre(Interval())
Jacobi(a,b)=Jacobi(a,b,Interval())

spacescompatible(a::Jacobi,b::Jacobi)=a.a==b.a && a.b==b.b

function canonicalspace(S::Jacobi)
    if isinteger(S.a) && isinteger(S.b)
        Jacobi(0.,0.,domain(S))
    elseif isinteger(S.a+0.5) && isinteger(S.b+0.5)
        Chebyshev()
    else
        error("There is no canonical space for Jacobi with a="*string(S.a)*" and b="*string(S.b))
    end
end

#####
# jacobirecA/B/C is from dlmf:
# p_{n+1} = (A_n x + B_n)p_n - C_n p_{n-1}
#####
jacobirecA(α,β,k)=k==0&&((α+β==0)||(α+β==-1))?.5*(α+β)+1:(2k+α+β+1)*(2k+α+β+2)/(2*(k+1)*(k+α+β+1))
jacobirecB(α,β,k)=k==0&&((α+β==0)||(α+β==-1))?.5*(β-α):(α^2-β^2)*(2k+α+β+1)/(2*(k+1)*(k+α+β+1)*(2k+α+β))
jacobirecC(α,β,k)=(k+α)*(k+β)*(2k+α+β+2)/((k+1)*(k+α+β+1)*(2k+α+β))

#####
# jacobirecA/B/C is from dlmf:
# x p_{n-1} =γ_n p_{n-2} + α_n p_{n-1} +  p_n β_n   
#####

jacobirecγ(α,β,k)=jacobirecC(α,β,k-1)/jacobirecA(α,β,k-1)
jacobirecα(α,β,k)=-jacobirecB(α,β,k-1)/jacobirecA(α,β,k-1)
jacobirecβ(α,β,k)=1/jacobirecA(α,β,k-1)

for (REC,JREC) in ((:recα,:jacobirecα),(:recβ,:jacobirecβ),(:recγ,:jacobirecγ))
    @eval $REC(::Type,sp::Jacobi,k)=$JREC(sp.a,sp.b,k)  #TODO: implement typing
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
            v=Array(Float64,n)
            v[1]=1.
            v[2]=.5*(α-β+(2+α+β)*x)
            
            for k=2:n-1
                v[k+1]=((x-jacobirecα(α,β,k))*v[k] - jacobirecγ(α,β,k)*v[k-1])/jacobirecβ(α,β,k)
            end
        end
        v[r+1]
    end
end
jacobip(n::Integer,α,β,v::Number)=jacobip(n:n,α,β,v)[1]
jacobip(n::Range1,α,β,v::Vector)=hcat(map(x->jacobip(n,α,β,x),v)...).'
jacobip(n::Integer,α,β,v::Vector)=map(x->jacobip(n,α,β,x),v)
jacobip(n,S::Jacobi,v)=jacobip(n,S.a,S.b,v)





include("jacobitransform.jl")
include("JacobiOperators.jl")



for op in (:(Base.ones),:(Base.zeros))
    @eval ($op){T<:Number}(::Type{T},S::Jacobi)=Fun(($op)(T,1),S)
    @eval ($op)(S::Jacobi)=Fun(($op)(1),S)    
end

function identity_fun(J::Jacobi)
    @assert domain(J)==Interval()
    Fun([(J.b-J.a)/(2+J.a+J.b),2.0/(2+J.a+J.b)],J)
end