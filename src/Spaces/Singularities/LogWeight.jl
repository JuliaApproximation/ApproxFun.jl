"""
    LogWeight(β,α,s::Space)

represents a function on `-1..1` weighted by `log((1+x)^β*(1-x)^α)`.
For other domains, the weight is inferred by mapping to `-1..1`.
"""
struct LogWeight{S,DD,RR,T} <: WeightSpace{S,DD,RR}
    β::T
    α::T
    space::S
end


LogWeight(β,α,space) =
    LogWeight{typeof(space),domaintype(space),rangetype(space),promote_type(typeof(β),typeof(α))}(β,α,space)

spacescompatible(A::LogWeight,B::LogWeight)=A.β==B.β && A.α == B.α && spacescompatible(A.space,B.space)
canonicalspace(A::LogWeight)=A

logweight(β,α,x)=log((1+x)^β*(1-x)^α)
weight(sp::LogWeight,x)=logweight(sp.β,sp.α,tocanonical(sp,x))


setdomain(sp::LogWeight,d::Domain)=LogWeight(sp.β,sp.α,setdomain(sp.space,d))

function coefficients(f::AbstractVector,sp1::LogWeight,sp2::LogWeight)
    β,α=sp1.β,sp1.α
    c,d=sp2.β,sp2.α

    if isapprox(c,β) && isapprox(d,α)
        coefficients(f,sp1.space,sp2.space)
    else
        (Conversion(sp1,sp2)*f)
    end
end

for (OPrule,OP) in ((:maxspace_rule,:maxspace),(:union_rule,:union))
    @eval begin
        function $OPrule(A::LogWeight,B::LogWeight)
            if isapprox(A.β,B.β) && isapprox(A.α,B.α)
                LogWeight(A.β,A.α,$OP(A.space,B.space))
            else
                NoSpace()
            end
        end
        # there are no other comatible spaces
        # this is mostly to overcome typing bug
        # in spacestes
        $OPrule(::LogWeight,::JacobiWeight)=NoSpace()
        $OPrule(::JacobiWeight,::LogWeight)=NoSpace()
    end
end








# Same as JacobiWeight

# avoid redundency
function Multiplication(f::Fun{JacobiWeight{SS,DD,RR,TT},T}, S::LogWeight{LWS,DD,RR}) where {SS,LWS,DD<:IntervalDomain,RR,T,TT}
    M=Multiplication(Fun(space(f).space,f.coefficients),S)
    rsp=JacobiWeight(space(f).β,space(f).α,rangespace(M))
    MultiplicationWrapper(f,SpaceOperator(M,S,rsp))
end

function Multiplication(f::Fun{<:LogWeight},S::JacobiWeight{SS,DD}) where {SS,DD<:IntervalDomain}
    M=Multiplication(f,S.space)
    rsp=JacobiWeight(S.β,S.α,rangespace(M))
    MultiplicationWrapper(f,SpaceOperator(M,S,rsp))
end


function Multiplication(f::Fun,S::LogWeight)
    M=Multiplication(f,S.space)
    rsp=LogWeight(S.β,S.α,rangespace(M))
    MultiplicationWrapper(f,SpaceOperator(M,S,rsp))
end

function Multiplication(f::Fun{<:LogWeight},S::LogWeight)
    M=Multiplication(f,S.space)
    rsp=LogWeight(S.β,S.α,rangespace(M))
    MultiplicationWrapper(f,SpaceOperator(M,S,rsp))
end



function Multiplication(f::Fun{<:LogWeight},S::Space)
    M=Multiplication(Fun(space(f).space,f.coefficients),S)
    rsp=LogWeight(space(f).β,space(f).α,rangespace(M))
    MultiplicationWrapper(f,SpaceOperator(M,S,rsp))
end
