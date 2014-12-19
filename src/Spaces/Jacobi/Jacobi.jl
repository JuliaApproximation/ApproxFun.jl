export JacobiSpace,LegendreSpace

immutable JacobiSpace <: IntervalSpace
    a::Float64
    b::Float64
    domain::Union(IntervalDomain,AnyDomain)
end
LegendreSpace(domain)=JacobiSpace(0.,0.,domain)
LegendreSpace()=LegendreSpace(Interval())
JacobiSpace(a,b)=JacobiSpace(a,b,Interval())

spacescompatible(a::JacobiSpace,b::JacobiSpace)=a.a==b.a && a.b==b.b

function canonicalspace(S::JacobiSpace)
    if isinteger(S.a) && isinteger(S.b)
        JacobiSpace(0.,0.,domain(S))
    elseif isinteger(S.a+0.5) && isinteger(S.b+0.5)
        ChebyshevSpace()
    else
        error("There is no canonical space for Jacobi with a="*string(S.a)*" and b="*string(S.b))
    end
end


jacobirecA(α,β,k)=k==0&&((α+β==0)||(α+β==-1))?.5*(α+β)+1:(2k+α+β+1)*(2k+α+β+2)/(2*(k+1)*(k+α+β+1))
jacobirecB(α,β,k)=k==0&&((α+β==0)||(α+β==-1))?.5*(β-α):(α^2-β^2)*(2k+α+β+1)/(2*(k+1)*(k+α+β+1)*(2k+α+β))
jacobirecC(α,β,k)=(k+α)*(k+β)*(2k+α+β+2)/((k+1)*(k+α+β+1)*(2k+α+β))

jacobirecγ(α,β,k)=jacobirecC(α,β,k-1)/jacobirecA(α,β,k-1)
jacobirecα(α,β,k)=-jacobirecB(α,β,k-1)/jacobirecA(α,β,k-1)
jacobirecβ(α,β,k)=1/jacobirecA(α,β,k-1)


type JacobiRecurrenceOperator <: TridiagonalOperator{Float64}
    a::Float64
    b::Float64
end

function getdiagonalentry(J::JacobiRecurrenceOperator,k,j)
    if j==-1
        jacobirecγ(J.a,J.b,k)
    elseif j==0
        jacobirecα(J.a,J.b,k)
    else #j==1
        jacobirecβ(J.a,J.b,k)
    end
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
jacobip(n,S::JacobiSpace,v)=jacobip(n,S.a,S.b,v)





include("jacobitransform.jl")
include("JacobiOperators.jl")



for op in (:(Base.ones),:(Base.zeros))
    @eval ($op){T<:Number}(::Type{T},S::JacobiSpace)=Fun(($op)(T,1),S)
    @eval ($op)(S::JacobiSpace)=Fun(($op)(1),S)    
end

function identity_fun(J::JacobiSpace)
    @assert domain(J)==Interval()
    Fun([(J.b-J.a)/(2+J.a+J.b),2.0/(2+J.a+J.b)],J)
end