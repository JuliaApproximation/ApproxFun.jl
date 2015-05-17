


## Calculus

for (Func,Len) in ((:(Base.sum),:complexlength),(:linesum,:length))
    @eval begin
        function $Func(f::Fun{JacobiWeight{Chebyshev}})
            d,α,β,n=domain(f),f.space.α,f.space.β,length(f)
            if α ≤ -1.0 || β ≤ -1.0
                fs = Fun(f.coefficients,f.space.space)
                return Inf*0.5*$Len(d)*(sign(fs[d.a])+sign(fs[d.b]))/2
            elseif α == β == -0.5
                return 0.5*$Len(d)*f.coefficients[1]*π
            elseif α == β == 0.5
                return 0.5*$Len(d)*(n ≤ 2 ? f.coefficients[1]/2 : f.coefficients[1]/2 - f.coefficients[3]/4)*π
            elseif α == 0.5 && β == -0.5
                return 0.5*$Len(d)*(n == 1 ? f.coefficients[1] : f.coefficients[1] + f.coefficients[2]/2)*π
            elseif α == -0.5 && β == 0.5
                return 0.5*$Len(d)*(n == 1 ? f.coefficients[1] : f.coefficients[1] - f.coefficients[2]/2)*π
            else
                c = zeros(eltype(f),n)
                c[1] = 2.^(α+β+1)*gamma(α+1)*gamma(β+1)/gamma(α+β+2)
                if n > 1
                    c[2] = c[1]*(α-β)/(α+β+2)
                    for i=1:n-2
                        c[i+2] = (2(α-β)*c[i+1]-(α+β-i+2)*c[i])/(α+β+i+2)
                    end
                end
                return 0.5*$Len(d)*dotu(f.coefficients,c)
            end
        end
        $Func{PS<:PolynomialSpace}(f::Fun{JacobiWeight{PS}})=$Func(Fun(f,
                                                                       JacobiWeight(space(f).α,
                                                                                    space(f).β,
                                                                                    Chebyshev(domain(f)))))
    end
end


function differentiate{J<:JacobiWeight}(f::Fun{J})
    S=f.space
    d=domain(f)
    ff=Fun(f.coefficients,S.space)
    if S.α==S.β==0
        u=differentiate(ff)
        Fun(u.coefficients,JacobiWeight(0.,0.,space(u)))
    elseif S.α==0
        x=Fun(identity,d)
        M=tocanonical(d,x)
        Mp=tocanonicalD(d,d.a)
        u=-Mp*S.β*ff +(1-M).*differentiate(ff)
        Fun(u.coefficients,JacobiWeight(0.,S.β-1,space(u)))
    elseif S.β==0
        x=Fun(identity,d)
        M=tocanonical(d,x)
        Mp=tocanonicalD(d,d.a)
        u=Mp*S.α*ff +(1+M).*differentiate(ff)
        Fun(u.coefficients,JacobiWeight(S.α-1,0.,space(u)))
    else
        x=Fun(identity,d)
        M=tocanonical(d,x)
        Mp=tocanonicalD(d,d.a)
        u=(Mp*S.α)*(1-M).*ff- (Mp*S.β)*(1+M).*ff +(1-M.^2).*differentiate(ff)
        Fun(u.coefficients,JacobiWeight(S.α-1,S.β-1,space(u)))
    end
end



function integrate{J<:JacobiWeight}(f::Fun{J})
    S=space(f)
    # we integrate by solving u'=f
    D=Derivative(S)
    if S.α==0 || S.β==0
        D\f
    else
        s=sum(f)
        if isapprox(s,0.)
            D\f
        else
            # we normalized so it sums to zero, and so backslash works
            w=Fun(x->exp(-40x^2),81)
            w1=Fun(coefficients(w),S)
            w2=Fun(x->w1[x],domain(w1))
            c=s/sum(w1)
            v=f-w1*c
            (c*integrate(w2))⊕linsolve(D,v;tolerance=100eps())
        end
    end
end

function Base.cumsum{J<:JacobiWeight}(f::Fun{J})
    g=integrate(f)

    S=space(f)

    if S.α==0 && S.β==0
        g-first(g)
    elseif S.α>-1
        Fun(-first(g),domain(g))⊕g
    else
        warn("Function is not integrable at left endpoint.  Returning a non-normalized indefinite integral.")
        g
    end
end


## Operators


function Derivative(S::JacobiWeight)
    d=domain(S)

    if S.α==S.β==0
        DerivativeWrapper(SpaceOperator(Derivative(S.space),S,JacobiWeight(0.,0.,rangespace(Derivative(S.space)))),1)
    elseif S.α==0
        x=Fun(identity,d)
        M=tocanonical(d,x)
        Mp=isa(d,Interval)?tocanonicalD(d,d.a):Fun(tocanonicalD(d,x),S.space) #TODO hack for Ray, which returns JacobiWeight but doesn't need to
        DD=(-Mp*S.β) +(1-M)*Derivative(S.space)
        DerivativeWrapper(SpaceOperator(DD,S,JacobiWeight(0.,S.β-1,rangespace(DD))),1)
    elseif S.β==0
        x=Fun(identity,d)
        M=tocanonical(d,x)
        Mp=isa(d,Interval)?tocanonicalD(d,d.a):Fun(tocanonicalD(d,x),S.space) #TODO hack for Ray, which returns JacobiWeight but doesn't need to
        DD=(Mp*S.α) +(1+M)*Derivative(S.space)
        DerivativeWrapper(SpaceOperator(DD,S,JacobiWeight(S.α-1,0.,rangespace(DD))),1)
    else
        x=Fun(identity,d)
        M=tocanonical(d,x)
        Mp=isa(d,Interval)?tocanonicalD(d,d.a):Fun(tocanonicalD(d,x),S.space) #TODO hack for Ray, which returns JacobiWeight but doesn't need to
        DD=(Mp*S.α)*(1-M) - (Mp*S.β)*(1+M) +(1-M.^2)*Derivative(S.space)
        DerivativeWrapper(SpaceOperator(DD,S,JacobiWeight(S.α-1,S.β-1,rangespace(DD))),1)
    end

end

function Derivative(S::JacobiWeight,k::Integer)
    if k==1
        Derivative(S)
    else
        D=Derivative(S)
        DerivativeWrapper(TimesOperator(Derivative(rangespace(D),k-1).op,D.op),k)
    end
end




## Multiplication

#Left multiplication. Here, S is considered the domainspace and we determine rangespace accordingly.

function Multiplication{D<:JacobiWeight,T}(f::Fun{D,T},S::JacobiWeight)
    M=Multiplication(Fun(f.coefficients,space(f).space),S.space)
    rsp=canonicalspace(JacobiWeight(space(f).α+S.α,space(f).β+S.β,rangespace(M)))
    MultiplicationWrapper(f,SpaceOperator(M,S,rsp))
end

function Multiplication{D,T}(f::Fun{D,T},S::JacobiWeight)
    M=Multiplication(f,S.space)
    rsp=JacobiWeight(S.α,S.β,rangespace(M))
    MultiplicationWrapper(f,SpaceOperator(M,S,rsp))
end

function Multiplication{D<:JacobiWeight,T}(f::Fun{D,T},S::IntervalSpace)
    M=Multiplication(Fun(f.coefficients,space(f).space),S)
    rsp=JacobiWeight(space(f).α,space(f).β,rangespace(M))
    MultiplicationWrapper(f,SpaceOperator(M,S,rsp))
end

#Right multiplication. Here, S is considered the rangespace and we determine domainspace accordingly.

function Multiplication{D<:JacobiWeight,T}(S::JacobiWeight,f::Fun{D,T})
    M=Multiplication(Fun(f.coefficients,space(f).space),S.space)
    dsp=canonicalspace(JacobiWeight(S.α-space(f).α,S.β-space(f).β,rangespace(M)))
    MultiplicationWrapper(f,SpaceOperator(M,dsp,S))
end

function Multiplication{D,T}(S::JacobiWeight,f::Fun{D,T})
    M=Multiplication(f,S.space)
    dsp=JacobiWeight(S.α,S.β,rangespace(M))
    MultiplicationWrapper(f,SpaceOperator(M,dsp,S))
end

function Multiplication{D<:JacobiWeight,T}(S::IntervalSpace,f::Fun{D,T})
    M=Multiplication(Fun(f.coefficients,space(f).space),S)
    dsp=JacobiWeight(-space(f).α,-space(f).β,rangespace(M))
    MultiplicationWrapper(f,SpaceOperator(M,dsp,S))
end

## Conversion

maxspace(A::JacobiWeight,B::JacobiWeight)=JacobiWeight(min(A.α,B.α),min(A.β,B.β),maxspace(A.space,B.space))
maxspace(A::IntervalSpace,B::JacobiWeight)=maxspace(JacobiWeight(0.,0.,A),B)
maxspace(A::JacobiWeight,B::IntervalSpace)=maxspace(A,JacobiWeight(0.,0.,B))

isapproxinteger(x)=isapprox(x,round(Int,x))

# return the space that has banded Conversion to the other, or NoSpace
conversion_rule{n,S<:FunctionSpace,IS<:IntervalSpace}(A::SliceSpace{n,1,S,RealBasis},B::JacobiWeight{IS})=error("Not implemented")
conversion_rule(A::JacobiWeight,B::JacobiWeight)=JacobiWeight(max(A.α,B.α),max(A.β,B.β),conversion_type(A.space,B.space))
#conversion_rule(A::IntervalSpace,B::JacobiWeight)=conversion_type(JacobiWeight(0,0,A),B)
conversion_rule(A::JacobiWeight,B::IntervalSpace)=conversion_type(A,JacobiWeight(0,0,B))



function Conversion(A::JacobiWeight,B::JacobiWeight)
    @assert isapproxinteger(A.α-B.α) && isapproxinteger(A.β-B.β)

    if isapprox(A.α,B.α) && isapprox(A.β,B.β)
        SpaceOperator(Conversion(A.space,B.space),A,B)
    elseif A.space==B.space
        @assert A.α≥B.α&&A.β≥B.β
        d=domain(A)
        x=Fun(identity,d)
        M=tocanonical(d,x)
        m=(1+M).^round(Int,A.α-B.α).*(1-M).^round(Int,A.β-B.β)
        MC=Multiplication(m,B.space)
        # The following is just a safety check
        @assert rangespace(MC) == B.space
        SpaceOperator(MC,A,B)# Wrap the operator with the correct spaces
    else
        @assert A.α≥B.α&&A.β≥B.β
        d=domain(A)
        x=Fun(identity,d)
        M=tocanonical(d,x)
        C=Conversion(A.space,B.space)
        m=(1+M).^round(Int,A.α-B.α).*(1-M).^round(Int,A.β-B.β)
        MC=TimesOperator(Multiplication(m,B.space),C)
        # The following is just a safety check
        @assert rangespace(MC) == B.space
        SpaceOperator(MC,A,B)
    end
end

Conversion(A::IntervalSpace,B::JacobiWeight)=ConversionWrapper(
    SpaceOperator(
        Conversion(JacobiWeight(0,0,A),B),
        A,B))
Conversion(A::JacobiWeight,B::IntervalSpace)=ConversionWrapper(
    SpaceOperator(
        Conversion(A,JacobiWeight(0,0,B)),
        A,B))






## Evaluation

function  Base.getindex{J<:JacobiWeight}(op::Evaluation{J,Bool},kr::Range)
    S=op.space
    @assert op.order<=1
    d=domain(op)

    if op.x
        @assert S.β>=0
        if S.β==0
            if op.order==0
                2^S.α*getindex(Evaluation(S.space,op.x),kr)
            else #op.order ===1
                @assert isa(d,Interval)
                2^S.α*getindex(Evaluation(S.space,op.x,1),kr)+(tocanonicalD(d,d.a)*S.α*2^(S.α-1))*getindex(Evaluation(S.space,op.x),kr)
            end
        else
            @assert op.order==0
            zeros(kr)
        end
    else
        @assert S.α>=0
        if S.α==0
            if op.order==0
                2^S.β*getindex(Evaluation(S.space,op.x),kr)
            else #op.order ===1
                @assert isa(d,Interval)
                2^S.β*getindex(Evaluation(S.space,op.x,1),kr)-(tocanonicalD(d,d.a)*S.β*2^(S.β-1))*getindex(Evaluation(S.space,op.x),kr)
            end
        else
            @assert op.order==0
            zeros(kr)
        end
    end
end


## Definite Integral

for (Func,Len) in ((:DefiniteIntegral,:complexlength),(:DefiniteLineIntegral,:length))
    @eval begin

        function getindex{λ,T}(Σ::$Func{JacobiWeight{Ultraspherical{λ}},T},kr::Range)
            dsp = domainspace(Σ)
            d = domain(Σ)
            @assert isa(d,Interval)
            @assert dsp.α==dsp.β==λ-0.5

            C = $Len(d)/2

            promote_type(T,typeof(C))[k == 1? C*gamma(λ+one(T)/2)*gamma(one(T)/2)/gamma(λ+one(T)) : zero(T) for k=kr]
        end

        datalength{λ}(Σ::$Func{JacobiWeight{Ultraspherical{λ}}})=1
    end
end
