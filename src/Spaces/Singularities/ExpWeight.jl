export ExpWeight

immutable ExpWeight{S,FF,DD} <: WeightSpace{S,RealBasis,DD,1}
    exponent::FF
    space::S
end

function ExpWeight(q,sp)
    @assert domain(q) == domain(sp)
    ExpWeight{typeof(sp),typeof(q),typeof(domain(q))}(q,sp)
end

weight(sp::ExpWeight,x) = exp(sp.exponent(x))
spacescompatible(a::ExpWeight,b::ExpWeight) = a.exponent ≈ b.exponent && spacescompatible(a.space,b.space)
setdomain(sp::ExpWeight,d::Domain) = ExpWeight(setdomain(sp.exponent,d),setdomain(sp.space,d))

for (RULE,TYP) in ((:conversion_rule,:conversion_type),(:maxspace_rule,:maxspace))
    @eval function $RULE(a::ExpWeight,b::ExpWeight)
        if a.exponent ≈ b.exponent
            sp=$TYP(a.space,b.space)
            if !isa(sp,NoSpace)
                return ExpWeight(a.exponent,sp)
            end
        end
        NoSpace()
    end
end

function Derivative(wsp::ExpWeight)
    q=wsp.exponent
    sp=wsp.space
    D=Derivative(sp) + q'
    DerivativeWrapper(SpaceOperator(D,wsp,ExpWeight(q,rangespace(D))),1)
end

function Derivative(S::ExpWeight,k::Integer)
    if k==1
        Derivative(S)
    else
        D=Derivative(S)
        DerivativeWrapper(TimesOperator(Derivative(rangespace(D),k-1),D),k)
    end
end

function Multiplication(f::Fun,sp::ExpWeight)
    M=Multiplication(f,sp.space)
    MultiplicationWrapper(f,SpaceOperator(M,sp,ExpWeight(sp.exponent,rangespace(M))))
end

function Conversion(a::ExpWeight,b::ExpWeight)
    @assert a.exponent ≈ b.exponent
    ConversionWrapper(SpaceOperator(Conversion(a.space,b.space),a,b))
end
