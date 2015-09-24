"""
 LogWeight
 represents a function on [-1,1] weighted by log((1+x)^α*(1-x)^β)
"""
immutable LogWeight{S,DD} <: WeightSpace{RealBasis,DD,1}
    α::Float64
    β::Float64
    space::S
end
LogWeight(α,β,space)=LogWeight{typeof(space),typeof(domain(space))}(α,β,space)

spacescompatible(A::LogWeight,B::LogWeight)=A.α==B.α && A.β == B.β && spacescompatible(A.space,B.space)
canonicalspace(A::LogWeight)=A

logweight(α,β,x)=log((1+x).^α.*(1-x).^β)
weight(sp::LogWeight,x)=logweight(sp.α,sp.β,tocanonical(sp,x))


setdomain(sp::LogWeight,d::Domain)=LogWeight(sp.α,sp.β,setdomain(sp.space,d))

function coefficients(f::Vector,sp1::LogWeight,sp2::LogWeight)
    α,β=sp1.α,sp1.β
    c,d=sp2.α,sp2.β

    if isapprox(c,α) && isapprox(d,β)
        coefficients(f,sp1.space,sp2.space)
    else
        (Conversion(sp1,sp2)*f)
    end
end

function maxspace_rule(A::LogWeight,B::LogWeight)
    if isapprox(A.α,B.α) && isapprox(A.β,B.β)
        LogWeight(A.α,A.β,maxspace(A.space,B.space))
    else
        NoSpace()
    end
end

# there are no other comatible spaces
# this is mostly to overcome typing bug
# in spacestes
maxspace_rule(::LogWeight,::JacobiWeight)=NoSpace()
maxspace_rule(::JacobiWeight,::LogWeight)=NoSpace()






# Same as JacobiWeight

# avoid redundency
function Multiplication{SS,LWS,DD<:Interval,T}(f::Fun{JacobiWeight{SS,DD},T},S::LogWeight{LWS,DD})
    M=Multiplication(Fun(f.coefficients,space(f).space),S)
    rsp=JacobiWeight(space(f).α,space(f).β,rangespace(M))
    MultiplicationWrapper(f,SpaceOperator(M,S,rsp))
end


function Multiplication(f::Fun,S::LogWeight)
    M=Multiplication(f,S.space)
    rsp=LogWeight(S.α,S.β,rangespace(M))
    MultiplicationWrapper(f,SpaceOperator(M,S,rsp))
end
