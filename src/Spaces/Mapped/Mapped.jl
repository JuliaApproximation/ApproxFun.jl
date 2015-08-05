export MappedSpace

##Mapped spaces

#Typing D as Domain was causing issues

type MappedSpace{S<:FunctionSpace,D,T} <: FunctionSpace{T,1}
    domain::D
    space::S
    MappedSpace(d::D,sp::S)=new(d,sp)
    MappedSpace(d::D)=new(d,S(canonicaldomain(d)))
    MappedSpace()=new(D(),S())
end


spacescompatible(a::MappedSpace,b::MappedSpace)=spacescompatible(a.space,b.space)&&domainscompatible(a,b)

MappedSpace{D<:Domain,T}(d::D,s::FunctionSpace{T})=MappedSpace{typeof(s),D,T}(d,s)




Space(d::Domain)=MappedSpace(d,Space(canonicaldomain(d)))


domain(S::MappedSpace)=S.domain
canonicalspace(S::MappedSpace)=MappedSpace(S.domain,canonicalspace(S.space))

points(d::MappedSpace,n)=fromcanonical(d,points(d.space,n))

## Construction

Base.ones{T<:Number}(::Type{T},S::MappedSpace)=Fun(ones(T,S.space).coefficients,S)
plan_transform(S::MappedSpace,vals::Vector)=plan_transform(S.space,vals)
plan_itransform(S::MappedSpace,cfs::Vector)=plan_itransform(S.space,cfs)
transform(S::MappedSpace,vals::Vector,plan...)=transform(S.space,vals,plan...)
itransform(S::MappedSpace,cfs::Vector,plan...)=itransform(S.space,cfs,plan...)
evaluate{SS,DD,T,TT}(f::Fun{MappedSpace{SS,DD,TT},T},x)=evaluate(Fun(coefficients(f),space(f).space),mappoint(domain(f),domain(space(f).space),x))


for op in (:(Base.first),:(Base.last))
    @eval $op{S<:MappedSpace}(f::Fun{S})=$op(Fun(coefficients(f),space(f).space))
end


## Integration


function integrate{LS,RR,T,TT}(f::Fun{MappedSpace{LS,RR,TT},T})
    fc=Fun(f.coefficients,space(f).space)
    x=Fun(identity,domain(fc))
    Mp=fromcanonicalD(f,x)
    g=integrate(fc*Mp)
    Fun(g.coefficients,MappedSpace(domain(f),space(g)))
end

function Base.sum{LS,RR,T,TT}(f::Fun{MappedSpace{LS,RR,TT},T})
    fc=Fun(f.coefficients,space(f).space)
    x=Fun(identity,domain(fc))
    Mp=fromcanonicalD(f,x)
    sum(fc*Mp)
end


function linesum{LS,RR,T,TT}(f::Fun{MappedSpace{LS,RR,TT},T})
    fc=Fun(f.coefficients,space(f).space)
    x=Fun(identity,domain(fc))
    Mp=fromcanonicalD(f,x)
    linesum(fc*abs(Mp))
end



## identity

function identity_fun{SS,DD,DDT}(S::MappedSpace{SS,DD,DDT})
    sf=fromcanonical(S,Fun(identity,domain(S.space)))
    Fun(coefficients(sf),MappedSpace(S.domain,space(sf)))
end

union_rule(A::ConstantSpace,B::MappedSpace)=MappedSpace(domain(B),union(A,B.space))

## Operators

function Evaluation(S1::MappedSpace,x::Bool,order::Integer)
    @assert order==0
    EvaluationWrapper(S1,x,order,Evaluation(S1.space,x,order))
end

Conversion(S1::MappedSpace,S2::MappedSpace)=ConversionWrapper(
    SpaceOperator(Conversion(S1.space,S2.space),
        S1,S2))

Conversion(S1::ConstantSpace,S2::MappedSpace)=ConversionWrapper(
    SpaceOperator(Conversion(S1,S2.space),
        S1,S2))

# Conversion is induced from canonical space
for OP in (:conversion_rule,:maxspace)
    @eval begin
        function $OP(S1::MappedSpace,S2::MappedSpace)
            @assert domain(S1)==domain(S2)
            cr=$OP(S1.space,S2.space)
            MappedSpace(domain(S1),cr)
        end
        function $OP(S1::ConstantSpace,S2::MappedSpace)
            cr=$OP(S1,S2.space)
            if isa(cr,ConstantSpace)
                cr
            else
                MappedSpace(domain(S2),cr)
            end
        end
    end
end

# Multiplication is the same as unmapped space
function Multiplication{MS<:MappedSpace,T}(f::Fun{MS,T},S::MappedSpace)
    d=domain(f)
    @assert d==domain(S)
    mf=Fun(coefficients(f),space(f).space)  # project f
    M=Multiplication(mf,S.space)
    MultiplicationWrapper(f,SpaceOperator(M,
        MappedSpace(d,domainspace(M)),
        MappedSpace(d,rangespace(M))
    ))
end


# Use tocanonicalD to find the correct derivative
function Derivative(S::MappedSpace,order::Int)
    x=Fun(identity,S)
    D1=Derivative(S.space)
    DS=SpaceOperator(D1,S,MappedSpace(domain(S),rangespace(D1)))
    M=Multiplication(Fun(tocanonicalD(S,x),S),DS|>rangespace)
    D=DerivativeWrapper(M*DS,1)
    if order==1
        D
    else
        Derivative(rangespace(D),order-1)*D
    end
end


function Integral(sp::MappedSpace,k::Integer)
    if k > 1
        Q=Integral(sp,1)
        IntegralWrapper(TimesOperator(Integral(rangespace(Q),k-1),Q),k)
    else # k==1
        csp=sp.space

        x=Fun(identity,csp)
        M=Multiplication(fromcanonicalD(sp,x),csp)
        Q=Integral(rangespace(M))*M
        IntegralWrapper(SpaceOperator(Q,sp,MappedSpace(sp.domain,rangespace(Q))),1)
    end
end


function DefiniteIntegral(sp::MappedSpace)
    x=Fun(domain(sp.space))
    M=Multiplication(fromcanonicalD(sp,x),sp.space)
    DefiniteIntegralWrapper(SpaceFunctional(DefiniteIntegral(rangespace(M))*M,sp))
end

function DefiniteLineIntegral(sp::MappedSpace)
    x=Fun(domain(sp.space))
    M=Multiplication(abs(fromcanonicalD(sp,x)),sp.space)
    DefiniteLineIntegralWrapper(SpaceFunctional(DefiniteIntegral(rangespace(M))*M,sp))
end


include("mappedinfdomains.jl")
include("CurveSpace.jl")
