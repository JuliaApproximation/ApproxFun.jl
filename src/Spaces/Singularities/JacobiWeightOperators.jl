


## Calculus

Base.sum(f::Fun{JacobiWeight{C,DD,RR,TT}}) where {C<:Chebyshev,DD,RR,TT} =
    sum(setdomain(f,canonicaldomain(f))*fromcanonicalD(f))
linesum(f::Fun{JacobiWeight{C,DD,RR,TT}}) where {C<:Chebyshev,DD,RR,TT} =
    linesum(setdomain(f,canonicaldomain(f))*abs(fromcanonicalD(f)))

for (Func,Len) in ((:(Base.sum),:complexlength),(:linesum,:arclength))
    @eval begin
        function $Func(f::Fun{JacobiWeight{C,DD,RR,TT}}) where {C<:Chebyshev,DD<:Segment,RR,TT}
            tol=1e-10
            d,β,α,n=domain(f),f.space.β,f.space.α,ncoefficients(f)
            g=Fun(space(f).space,f.coefficients)
            if β ≤ -1.0 && abs(first(g))≤tol
                $Func(increase_jacobi_parameter(-1,f))
            elseif α ≤ -1.0 && abs(last(g))≤tol
                $Func(increase_jacobi_parameter(+1,f))
            elseif β ≤ -1.0 || α ≤ -1.0
                fs = Fun(f.space.space,f.coefficients)
                return Inf*0.5*$Len(d)*(sign(fs(d.a))+sign(fs(d.b)))/2
            elseif β == α == -0.5
                return 0.5*$Len(d)*f.coefficients[1]*π
            elseif β == α == 0.5
                return 0.5*$Len(d)*(n ≤ 2 ? f.coefficients[1]/2 : f.coefficients[1]/2 - f.coefficients[3]/4)*π
            elseif β == 0.5 && α == -0.5
                return 0.5*$Len(d)*(n == 1 ? f.coefficients[1] : f.coefficients[1] + f.coefficients[2]/2)*π
            elseif β == -0.5 && α == 0.5
                return 0.5*$Len(d)*(n == 1 ? f.coefficients[1] : f.coefficients[1] - f.coefficients[2]/2)*π
            else
                c = zeros(cfstype(f),n)
                c[1] = 2^(β+α+1)*gamma(β+1)*gamma(α+1)/gamma(β+α+2)
                if n > 1
                    c[2] = c[1]*(β-α)/(β+α+2)
                    for i=1:n-2
                        c[i+2] = (2(β-α)*c[i+1]-(β+α-i+2)*c[i])/(β+α+i+2)
                    end
                end
                return 0.5*$Len(d)*dotu(f.coefficients,c)
            end
        end
        $Func(f::Fun{JacobiWeight{PS,DD,RR,TT}}) where {PS<:PolynomialSpace,DD,RR,TT} =
            $Func(Fun(f,JacobiWeight(space(f).β, space(f).α, Chebyshev(domain(f)))))
    end
end


function differentiate(f::Fun{JacobiWeight{SS,DD,RR,TT}}) where {SS,DD<:Segment,RR,TT}
    S=f.space
    d=domain(f)
    ff=Fun(S.space,f.coefficients)
    if S.β==S.α==0
        u=differentiate(ff)
        Fun(JacobiWeight(0.,0.,space(u)),u.coefficients)
    elseif S.β==0
        x=Fun(identity,d)
        M=tocanonical(d,x)
        Mp=tocanonicalD(d,d.a)
        u=-Mp*S.α*ff +(1-M).*differentiate(ff)
        Fun(JacobiWeight(0.,S.α-1,space(u)),u.coefficients)
    elseif S.α==0
        x=Fun(identity,d)
        M=tocanonical(d,x)
        Mp=tocanonicalD(d,d.a)
        u=Mp*S.β*ff +(1+M).*differentiate(ff)
        Fun(JacobiWeight(S.β-1,0.,space(u)),u.coefficients)
    else
        x=Fun(identity,d)
        M=tocanonical(d,x)
        Mp=tocanonicalD(d,d.a)
        u=(Mp*S.β)*(1-M).*ff- (Mp*S.α)*(1+M).*ff +(1-M.^2).*differentiate(ff)
        Fun(JacobiWeight(S.β-1,S.α-1,space(u)),u.coefficients)
    end
end



function integrate(f::Fun{JacobiWeight{SS,DD,RR,TT}}) where {SS,DD<:Segment,RR,TT}
    S=space(f)
    # we integrate by solving u'=f
    tol=1e-10
    g=Fun(S.space,f.coefficients)
    if S.β ≈ 0 && S.α ≈ 0
        integrate(g)
    elseif S.β ≤ -1 && abs(first(g)) ≤ tol
        integrate(increase_jacobi_parameter(-1,f))
    elseif S.α ≤ -1 && abs(last(g)) ≤ tol
        integrate(increase_jacobi_parameter(+1,f))
    elseif S.β ≈ -1 && S.α ≈ -1
        error("Implement")
    elseif S.β ≈ -1 && S.α ≈ 0
        p=first(g)  # first value without weight
        fp = Fun(f-Fun(S,[p]),S.space)  # Subtract out right value and divide singularity via conversion
        d=domain(f)
        Mp=tocanonicalD(d,d.a)
        integrate(fp) ⊕ Fun(LogWeight(1.,0.,S.space),[p/Mp])
    elseif S.β ≈ -1 && S.α > 0 && isapproxinteger(S.α)
        # convert to zero case and integrate
        integrate(Fun(f,JacobiWeight(S.β,0.,S.space)))
    elseif S.α ≈ -1 && S.β ≈ 0
        p=last(g)  # last value without weight
        fp = Fun(f-Fun(S,[p]),S.space)  # Subtract out right value and divide singularity via conversion
        d=domain(f)
        Mp=tocanonicalD(d,d.a)
        integrate(fp) ⊕ Fun(LogWeight(zero(TT),one(TT),S.space),[-p/Mp])
    elseif isapprox(S.α,-1) && S.β > 0 && isapproxinteger(S.β)
        # convert to zero case and integrate
        integrate(Fun(f,JacobiWeight(zero(TT),S.α,S.space)))
    elseif S.β ≈ 0
        D = Derivative(JacobiWeight(S.β, S.α+1, S.space))
        D\f   # this happens to pick out a smooth solution
    elseif S.α ≈ 0
        D = Derivative(JacobiWeight(S.β+1, S.α, S.space))
        D\f   # this happens to pick out a smooth solution
    elseif isapproxinteger(S.β) || isapproxinteger(S.α)
        D = Derivative(JacobiWeight(S.β+1, S.α+1, S.space))
        D\f   # this happens to pick out a smooth solution
    else
        s=sum(f)
        if abs(s)<1E-14
            D=Derivative(JacobiWeight(S.β+1, S.α+1, S.space))
            \(D,f; tolerance=1E-14)  # if the sum is 0 we don't get step-like behaviour
        else
            # we normalized so it sums to zero, and so backslash works
            w = Fun(x->exp(-40x^2),81)
            w1 = Fun(S,coefficients(w))
            w2 = Fun(x->w1(x),domain(w1))
            c  = s/sum(w1)
            v  = f-w1*c
            (c*integrate(w2)) ⊕ integrate(v)
        end
    end
end

function Base.cumsum(f::Fun{JacobiWeight{SS,DD,RR,TT}}) where {SS,DD<:Segment,RR,TT}
    g=integrate(f)
    S=space(f)

    if (S.β==0 && S.α==0) || S.β>-1
        g-first(g)
    else
        @warn "Function is not integrable at left endpoint.  Returning a non-normalized indefinite integral."
        g
    end
end


## Operators


function jacobiweightDerivative(S::JacobiWeight{SS,DDD}) where {SS,DDD<:Segment}
    d=domain(S)

    if d!=Segment()
        # map to canonical
        Mp=fromcanonicalD(d,d.a)
        DD=jacobiweightDerivative(setdomain(S,Segment()))

        return DerivativeWrapper(SpaceOperator(DD.op.op,S,setdomain(rangespace(DD),d))/Mp,1)
    end


    if S.β==S.α==0
        DerivativeWrapper(SpaceOperator(Derivative(S.space),S,JacobiWeight(0.,0.,rangespace(Derivative(S.space)))),1)
    elseif S.β==0
        w=Fun(JacobiWeight(0,1,ConstantSpace(d)),[1.])

        DD=-S.α + w*Derivative(S.space)
        rs=S.α==1 ? rangespace(DD) : JacobiWeight(0.,S.α-1,rangespace(DD))
        DerivativeWrapper(SpaceOperator(DD,S,rs),1)
    elseif S.α==0
        w=Fun(JacobiWeight(1,0,ConstantSpace(d)),[1.])

        DD=S.β + w*Derivative(S.space)
        rs=S.β==1 ? rangespace(DD) : JacobiWeight(S.β-1,0.,rangespace(DD))
        DerivativeWrapper(SpaceOperator(DD,S,rs),1)
    else
        w=Fun(JacobiWeight(1,1,ConstantSpace(d)),[1.])
        x=Fun()

        DD=S.β*(1-x) - S.α*(1+x) + w*Derivative(S.space)
        rs=S.β==1&&S.α==1 ? rangespace(DD) : JacobiWeight(S.β-1,S.α-1,rangespace(DD))
        DerivativeWrapper(SpaceOperator(DD,S,rs),1)
    end
end

Derivative(S::JacobiWeight{SS,DDD}) where {SS,DDD<:Segment} = jacobiweightDerivative(S)

function Derivative(S::JacobiWeight{SS,DD}, k::Integer) where {SS,DD<:Segment}
    if k==1
        Derivative(S)
    else
        D=Derivative(S)
        DerivativeWrapper(TimesOperator(Derivative(rangespace(D),k-1),D),k)
    end
end




## Multiplication

#Left multiplication. Here, S is considered the domainspace and we determine rangespace accordingly.

function Multiplication(f::Fun{JacobiWeight{S1,DD,RR,TT},T},S::JacobiWeight{S2,DD,RR}) where {S1,S2,DD<:IntervalDomain,RR,T,TT}
    M=Multiplication(Fun(space(f).space,f.coefficients),S.space)
    if space(f).β+S.β==space(f).α+S.α==0
        rsp=rangespace(M)
    else
        rsp=JacobiWeight(space(f).β+S.β,space(f).α+S.α,rangespace(M))
    end
    MultiplicationWrapper(f,SpaceOperator(M,S,rsp))
end

function Multiplication(f::Fun{D,T}, S::JacobiWeight{SS,DD}) where {D,T,SS,DD<:IntervalDomain}
    M=Multiplication(f,S.space)
    rsp=JacobiWeight(S.β,S.α,rangespace(M))
    MultiplicationWrapper(f,SpaceOperator(M,S,rsp))
end

function Multiplication(f::Fun{JacobiWeight{SS,ID,RR,TT},T},S::PolynomialSpace{ID}) where {SS,T,ID<:IntervalDomain,RR,TT}
    M=Multiplication(Fun(space(f).space,f.coefficients),S)
    rsp=JacobiWeight(space(f).β,space(f).α,rangespace(M))
    MultiplicationWrapper(f,SpaceOperator(M,S,rsp))
end

#Right multiplication. Here, S is considered the rangespace and we determine domainspace accordingly.

function Multiplication(S::JacobiWeight{SS,DD,RR},f::Fun{JacobiWeight{S2,DD,RR,TT},T}) where {DD<:IntervalDomain,RR,SS,S2,T,TT}
    M=Multiplication(Fun(space(f).space,f.coefficients),S.space)
    dsp=canonicalspace(JacobiWeight(S.β-space(f).β,S.α-space(f).α,rangespace(M)))
    MultiplicationWrapper(f,SpaceOperator(M,dsp,S))
end

function Multiplication(S::JacobiWeight{SS,DD},f::Fun{D,T}) where {D,SS,DD<:IntervalDomain,T}
    M=Multiplication(f,S.space)
    dsp=JacobiWeight(S.β,S.α,rangespace(M))
    MultiplicationWrapper(f,SpaceOperator(M,dsp,S))
end

# function Multiplication{D<:JacobiWeight,T,V,ID<:IntervalDomain}(S::Space{V,D},f::Fun{ID,T})
#     M=Multiplication(Fun(f.coefficients,space(f).space),S)
#     dsp=JacobiWeight(-space(f).β,-space(f).α,rangespace(M))
#     MultiplicationWrapper(f,SpaceOperator(M,dsp,S))
# end

## Conversion
for (OPrule,OP) in ((:maxspace_rule,:maxspace),(:union_rule,:union))
    @eval begin
        function $OPrule(A::JacobiWeight,B::JacobiWeight)
            if domainscompatible(A,B) && isapproxinteger(A.β-B.β) && isapproxinteger(A.α-B.α)
                ms=$OP(A.space,B.space)
                if min(A.β,B.β)==0.0 && min(A.α,B.α) == 0.0
                    return ms
                else
                    return JacobiWeight(min(A.β,B.β),min(A.α,B.α),ms)
                end
            end
            NoSpace()
        end
        $OPrule(A::JacobiWeight,B::Space{D}) where {D<:IntervalDomain} = $OPrule(A,JacobiWeight(0.,0.,B))
    end
end


for FUNC in (:hasconversion,:isconvertible)
    @eval begin
        $FUNC(A::JacobiWeight{S1,D},B::JacobiWeight{S2,D}) where {S1,S2,D<:IntervalDomain} =
            isapproxinteger(A.β-B.β) &&
            isapproxinteger(A.α-B.α) && A.β ≥ B.β && A.α ≥ B.α && $FUNC(A.space,B.space)

        $FUNC(A::JacobiWeight{S,D},B::Space{D}) where {S,D<:IntervalDomain} =
            $FUNC(A,JacobiWeight(0.,0.,B))
        $FUNC(B::Space{D},A::JacobiWeight{S,D}) where {S,D<:IntervalDomain} =
            $FUNC(JacobiWeight(0.,0.,B),A)
    end
end


# return the space that has banded Conversion to the other, or NoSpace

function conversion_rule(A::JacobiWeight,B::JacobiWeight)
    if isapproxinteger(A.β-B.β) && isapproxinteger(A.α-B.α)
        ct=conversion_type(A.space,B.space)
        ct==NoSpace() ? NoSpace() : JacobiWeight(max(A.β,B.β),max(A.α,B.α),ct)
    else
        NoSpace()
    end
end

conversion_rule(A::JacobiWeight,B::Space{D}) where {D<:IntervalDomain} = conversion_type(A,JacobiWeight(0,0,B))


# override defaultConversion instead of Conversion to avoid ambiguity errors
function defaultConversion(A::JacobiWeight{JS1,DD},B::JacobiWeight{JS2,DD}) where {JS1,JS2,DD<:IntervalDomain}
    @assert isapproxinteger(A.β-B.β) && isapproxinteger(A.α-B.α)

    if isapprox(A.β,B.β) && isapprox(A.α,B.α)
        ConversionWrapper(SpaceOperator(Conversion(A.space,B.space),A,B))
    else
        @assert A.β≥B.β&&A.α≥B.α
        # first check if a multiplication by JacobiWeight times B.space is overloaded
        # this is currently designed for Jacobi multiplied by (1-x), etc.
        βdif,αdif=round(Int,A.β-B.β),round(Int,A.α-B.α)
        d=domain(A)
        M=Multiplication(jacobiweight(βdif,αdif,d),
                         A.space)

        if rangespace(M)==JacobiWeight(βdif,αdif,A.space)
            # M is the default, so we should use multiplication by polynomials instead
            x=Fun(identity,d)
            y=mobius(d,x)   # we use mobius instead of tocanonical so that it works for Funs
            m=(1+y).^βdif.*(1-y).^αdif
            MC=promoterangespace(Multiplication(m,A.space),B.space)

            ConversionWrapper(SpaceOperator(MC,A,B))# Wrap the operator with the correct spaces
        else
            ConversionWrapper(SpaceOperator(promoterangespace(M,B.space),
                                            A,B))
        end
    end
end

defaultConversion(A::Space{D,R},B::JacobiWeight{JS,D}) where {JS,D<:IntervalDomain,R<:Real}=ConversionWrapper(
    SpaceOperator(
        Conversion(JacobiWeight(0,0,A),B),
        A,B))
defaultConversion(A::JacobiWeight{JS,D},B::Space{D,R}) where {JS,D<:IntervalDomain,R<:Real}=ConversionWrapper(
    SpaceOperator(
        Conversion(A,JacobiWeight(0,0,B)),
        A,B))






## Evaluation

function  Base.getindex(op::ConcreteEvaluation{<:JacobiWeight,typeof(first)},kr::AbstractRange)
    S=op.space
    @assert op.order ≤ 1
    d=domain(op)

    @assert S.β ≥ 0
    if S.β==0
        if op.order==0
            2^S.α*getindex(Evaluation(S.space,op.x),kr)
        else #op.order ===1
            @assert isa(d,IntervalDomain)
            2^S.α*getindex(Evaluation(S.space,op.x,1),kr)-(tocanonicalD(d,d.a)*S.α*2^(S.α-1))*getindex(Evaluation(S.space,op.x),kr)
        end
    else
        @assert op.order==0
        zeros(eltype(op), length(kr))
    end
end

function  Base.getindex(op::ConcreteEvaluation{<:JacobiWeight,typeof(last)},kr::AbstractRange)
    S=op.space
    @assert op.order<=1
    d=domain(op)

    @assert S.α>=0
    if S.α==0
        if op.order==0
            2^S.β*getindex(Evaluation(S.space,op.x),kr)
        else #op.order ===1
            @assert isa(d,IntervalDomain)
            2^S.β*getindex(Evaluation(S.space,op.x,1),kr)+
                (tocanonicalD(d,d.a)*S.β*2^(S.β-1))*getindex(Evaluation(S.space,op.x),kr)
        end
    else
        @assert op.order==0
        zeros(eltype(op), length(kr))
    end
end


## Definite Integral

for (Func,Len,Sum) in ((:DefiniteIntegral,:complexlength,:sum),(:DefiniteLineIntegral,:arclength,:linesum))
    ConcFunc = Meta.parse("Concrete"*string(Func))

    @eval begin
        $Func(S::JacobiWeight{SS,D}) where {SS,D<:Segment} = $ConcFunc(S)

        getindex(Σ::$ConcFunc,k::Integer) = eltype(Σ)($Sum(Fun(domainspace(Σ),[zeros(eltype(Σ),k-1);1])))

        function getindex(Σ::$ConcFunc{JacobiWeight{Ultraspherical{LT,D,R},D,R,TT},T},k::Integer) where {LT,D<:Segment,R,T,TT}
            λ = order(domainspace(Σ).space)
            dsp = domainspace(Σ)
            d = domain(Σ)
            C = $Len(d)/2

            if dsp.β==dsp.α==λ-0.5
                k == 1 ? convert(T,C*gamma(λ+one(T)/2)*gamma(one(T)/2)/gamma(λ+one(T))) : zero(T)
            else
                convert(T,$Sum(Fun(dsp,[zeros(T,k-1);1])))
            end
        end

        function getindex(Σ::$ConcFunc{JacobiWeight{Ultraspherical{LT,D,R},D,R,TT},T},kr::AbstractRange) where {LT,D<:Segment,R,T,TT}
            λ = order(domainspace(Σ).space)
            dsp = domainspace(Σ)
            d = domain(Σ)
            C = $Len(d)/2

            if dsp.β==dsp.α==λ-0.5
                T[k == 1 ? C*gamma(λ+one(T)/2)*gamma(one(T)/2)/gamma(λ+one(T)) : zero(T) for k=kr]
            else
                T[$Sum(Fun(dsp,[zeros(T,k-1);1])) for k=kr]
            end
        end

        function bandwidths(Σ::$ConcFunc{JacobiWeight{Ultraspherical{LT,D,R},D,R,TT}}) where {LT,D<:Segment,R,TT}
            λ = order(domainspace(Σ).space)
            β,α = domainspace(Σ).β,domainspace(Σ).α
            if β==α && isapproxinteger(β-0.5-λ) && λ ≤ ceil(Int,β)
                0,2*(ceil(Int,β)-λ)
            else
                0,∞
            end
        end


        function getindex(Σ::$ConcFunc{JacobiWeight{Chebyshev{D,R},D,R,TT},T},k::Integer) where {D<:Segment,R,T,TT}
            dsp = domainspace(Σ)
            d = domain(Σ)
            C = $Len(d)/2

            if dsp.β==dsp.α==-0.5
                k == 1 ? convert(T,C*π) : zero(T)
            else
                convert(T,$Sum(Fun(dsp,[zeros(T,k-1);1])))
            end
        end

        function getindex(Σ::$ConcFunc{JacobiWeight{Chebyshev{D,R},D,R,TT},T},kr::AbstractRange) where {D<:Segment,R,T,TT}
            dsp = domainspace(Σ)
            d = domain(Σ)
            C = $Len(d)/2

            if dsp.β==dsp.α==-0.5
                T[k == 1 ? C*π : zero(T) for k=kr]
            else
                T[$Sum(Fun(dsp,[zeros(T,k-1);1])) for k=kr]
            end
        end

        function bandwidths(Σ::$ConcFunc{JacobiWeight{Chebyshev{D,R},D,R,TT}}) where {D<:Segment,R,TT}
            β,α = domainspace(Σ).β,domainspace(Σ).α
            if β==α && isapproxinteger(β-0.5) && 0 ≤ ceil(Int,β)
                0,2ceil(Int,β)
            else
                0,∞
            end
        end
    end
end
