


## Calculus

Base.sum{C<:Chebyshev,DD}(f::Fun{JacobiWeight{C,DD}})=sum(setdomain(f,canonicaldomain(f))*fromcanonicalD(f))
linesum{C<:Chebyshev,DD}(f::Fun{JacobiWeight{C,DD}})=linesum(setdomain(f,canonicaldomain(f))*abs(fromcanonicalD(f)))

for (Func,Len) in ((:(Base.sum),:complexlength),(:linesum,:length))
    @eval begin
        function $Func{C<:Chebyshev,DD<:Interval}(f::Fun{JacobiWeight{C,DD}})
            tol=1e-10
            d,α,β,n=domain(f),f.space.α,f.space.β,length(f)
            g=Fun(f.coefficients,space(f).space)
            if α ≤ -1.0 && abs(first(g))≤tol
                $Func(increase_jacobi_parameter(-1,f))
            elseif β ≤ -1.0 && abs(last(g))≤tol
                $Func(increase_jacobi_parameter(+1,f))
            elseif α ≤ -1.0 || β ≤ -1.0
                fs = Fun(f.coefficients,f.space.space)
                return Inf*0.5*$Len(d)*(sign(fs(d.a))+sign(fs(d.b)))/2
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
        $Func{PS<:PolynomialSpace,DD}(f::Fun{JacobiWeight{PS,DD}})=$Func(Fun(f,
                                                                       JacobiWeight(space(f).α,
                                                                                    space(f).β,
                                                                                    Chebyshev(domain(f)))))
    end
end


function differentiate{J<:JacobiWeight,DD<:Interval}(f::Fun{J,DD})
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



function integrate{SS,DD<:Interval}(f::Fun{JacobiWeight{SS,DD}})
    S=space(f)
    # we integrate by solving u'=f
    D=Derivative(S)
    tol=1e-10
    g=Fun(f.coefficients,S.space)
    if isapprox(S.α,0.) && isapprox(S.β,0.)
        integrate(g)
    elseif S.α ≤ -1.0 && abs(first(g))≤tol
        integrate(increase_jacobi_parameter(-1,f))
    elseif S.β ≤ -1.0 && abs(last(g))≤tol
        integrate(increase_jacobi_parameter(+1,f))
    elseif isapprox(S.α,-1) && isapprox(S.β,-1)
        error("Implement")
    elseif isapprox(S.α,-1) && isapprox(S.β,0)
        p=first(g)  # first value without weight
        fp = Fun(f-Fun([p],S),S.space)  # Subtract out right value and divide singularity via conversion
        d=domain(f)
        Mp=tocanonicalD(d,d.a)
        integrate(fp)⊕Fun([p/Mp],LogWeight(1.,0.,S.space))
    elseif isapprox(S.α,-1) && S.β > 0 && isapproxinteger(S.β)
        # convert to zero case and integrate
        integrate(Fun(f,JacobiWeight(S.α,0.,S.space)))
    elseif isapprox(S.β,-1) && isapprox(S.α,0.)
        p=last(g)  # last value without weight
        fp = Fun(f-Fun([p],S),S.space)  # Subtract out right value and divide singularity via conversion
        d=domain(f)
        Mp=tocanonicalD(d,d.a)
        integrate(fp)⊕Fun([-p/Mp],LogWeight(0.,1.,S.space))
    elseif isapprox(S.β,-1) && S.α > 0 && isapproxinteger(S.α)
        # convert to zero case and integrate
        integrate(Fun(f,JacobiWeight(0.,S.β,S.space)))
    elseif isapprox(S.α,0) || isapprox(S.β,0)
        D\f   # this happens to pick out a smooth solution
    else
        s=sum(f)
        if abs(s)<1E-14
            linsolve(D,f;tolerance=1E-14)  # if the sum is 0 we don't get step-like behaviour
        else
            # we normalized so it sums to zero, and so backslash works
            w=Fun(x->exp(-40x^2),81)
            w1=Fun(coefficients(w),S)
            w2=Fun(x->w1(x),domain(w1))
            c=s/sum(w1)
            v=f-w1*c
            (c*integrate(w2))⊕linsolve(D,v;tolerance=100eps())
        end
    end
end

function Base.cumsum{SS,DD<:Interval}(f::Fun{JacobiWeight{SS,DD}})
    g=integrate(f)
    S=space(f)

    if (S.α==0 && S.β==0) || S.α>-1
        g-first(g)
    else
        warn("Function is not integrable at left endpoint.  Returning a non-normalized indefinite integral.")
        g
    end
end


## Operators


function Derivative{SS,DDD<:Interval}(S::JacobiWeight{SS,DDD})
    d=domain(S)

    if S.α==S.β==0
        DerivativeWrapper(SpaceOperator(Derivative(S.space),S,JacobiWeight(0.,0.,rangespace(Derivative(S.space)))),1)
    elseif S.α==0
        x=Fun(identity,d)
        M=tocanonical(d,x)
        Mp=isa(d,Interval)?tocanonicalD(d,d.a):Fun(tocanonicalD(d,x),S.space) #TODO hack for Ray, which returns JacobiWeight but doesn't need to
        DD=(-Mp*S.β) +(1-M)*Derivative(S.space)
        rs=S.β==1?rangespace(DD):JacobiWeight(0.,S.β-1,rangespace(DD))
        DerivativeWrapper(SpaceOperator(DD,S,rs),1)
    elseif S.β==0
        x=Fun(identity,d)
        M=tocanonical(d,x)
        Mp=isa(d,Interval)?tocanonicalD(d,d.a):Fun(tocanonicalD(d,x),S.space) #TODO hack for Ray, which returns JacobiWeight but doesn't need to
        DD=(Mp*S.α) +(1+M)*Derivative(S.space)
        rs=S.α==1?rangespace(DD):JacobiWeight(S.α-1,0.,rangespace(DD))
        DerivativeWrapper(SpaceOperator(DD,S,rs),1)
    else
        x=Fun(identity,d)
        M=tocanonical(d,x)
        Mp=isa(d,Interval)?tocanonicalD(d,d.a):Fun(tocanonicalD(d,x),S.space) #TODO hack for Ray, which returns JacobiWeight but doesn't need to
        DD=(Mp*S.α)*(1-M) - (Mp*S.β)*(1+M) +(1-M.^2)*Derivative(S.space)
        rs=S.α==1&&s.β==1?rangespace(DD):JacobiWeight(S.α-1,S.β-1,rangespace(DD))
        DerivativeWrapper(SpaceOperator(DD,S,rs),1)
    end

end

function Derivative{SS,DD<:Interval}(S::JacobiWeight{SS,DD},k::Integer)
    if k==1
        Derivative(S)
    else
        D=Derivative(S)
        DerivativeWrapper(TimesOperator(Derivative(rangespace(D),k-1),D),k)
    end
end




## Multiplication

#Left multiplication. Here, S is considered the domainspace and we determine rangespace accordingly.

function Multiplication{S1,S2,DD<:IntervalDomain,T}(f::Fun{JacobiWeight{S1,DD},T},S::JacobiWeight{S2,DD})
    M=Multiplication(Fun(f.coefficients,space(f).space),S.space)
    rsp=canonicalspace(JacobiWeight(space(f).α+S.α,space(f).β+S.β,rangespace(M)))
    MultiplicationWrapper(f,SpaceOperator(M,S,rsp))
end

function Multiplication{D,T,SS,DD<:IntervalDomain}(f::Fun{D,T},S::JacobiWeight{SS,DD})
    M=Multiplication(f,S.space)
    rsp=JacobiWeight(S.α,S.β,rangespace(M))
    MultiplicationWrapper(f,SpaceOperator(M,S,rsp))
end

function Multiplication{SS,T,V,ID<:IntervalDomain}(f::Fun{JacobiWeight{SS,ID},T},S::Space{V,ID})
    M=Multiplication(Fun(f.coefficients,space(f).space),S)
    rsp=JacobiWeight(space(f).α,space(f).β,rangespace(M))
    MultiplicationWrapper(f,SpaceOperator(M,S,rsp))
end

#Right multiplication. Here, S is considered the rangespace and we determine domainspace accordingly.

function Multiplication{DD<:IntervalDomain,SS,S2,T}(S::JacobiWeight{SS,DD},f::Fun{JacobiWeight{S2,DD},T})
    M=Multiplication(Fun(f.coefficients,space(f).space),S.space)
    dsp=canonicalspace(JacobiWeight(S.α-space(f).α,S.β-space(f).β,rangespace(M)))
    MultiplicationWrapper(f,SpaceOperator(M,dsp,S))
end

function Multiplication{D,SS,DD<:IntervalDomain,T}(S::JacobiWeight{SS,DD},f::Fun{D,T})
    M=Multiplication(f,S.space)
    dsp=JacobiWeight(S.α,S.β,rangespace(M))
    MultiplicationWrapper(f,SpaceOperator(M,dsp,S))
end

# function Multiplication{D<:JacobiWeight,T,V,ID<:IntervalDomain}(S::Space{V,D},f::Fun{ID,T})
#     M=Multiplication(Fun(f.coefficients,space(f).space),S)
#     dsp=JacobiWeight(-space(f).α,-space(f).β,rangespace(M))
#     MultiplicationWrapper(f,SpaceOperator(M,dsp,S))
# end

## Conversion

isapproxinteger(x)=isapprox(x,round(Int,x))

function maxspace_rule(A::JacobiWeight,B::JacobiWeight)
    if isapproxinteger(A.α-B.α) && isapproxinteger(A.β-B.β)
        ms=maxspace(A.space,B.space)
        if min(A.α,B.α)==0.&&min(A.β,B.β)==0.
            return ms
        else
            return JacobiWeight(min(A.α,B.α),min(A.β,B.β),ms)
        end
    end
    NoSpace()
end
maxspace_rule{T,D<:IntervalDomain}(A::JacobiWeight,B::Space{T,D})=maxspace(A,JacobiWeight(0.,0.,B))


hasconversion{S1,S2,D<:IntervalDomain}(A::JacobiWeight{S1,D},B::JacobiWeight{S2,D})=isapproxinteger(A.α-B.α) &&
    isapproxinteger(A.β-B.β) && A.α ≥ B.α && A.β ≥ B.β && hasconversion(A.space,B.space)


hasconversion{T,S,D<:IntervalDomain}(A::JacobiWeight{S,D},B::Space{T,D})=hasconversion(A,JacobiWeight(0.,0.,B))
hasconversion{T,S,D<:IntervalDomain}(B::Space{T,D},A::JacobiWeight{S,D})=hasconversion(JacobiWeight(0.,0.,B),A)



# return the space that has banded Conversion to the other, or NoSpace
conversion_rule{n,S<:Space}(A::SliceSpace{n,1,S,RealBasis},B::JacobiWeight)=error("Not implemented")
function conversion_rule(A::JacobiWeight,B::JacobiWeight)
    if isapproxinteger(A.α-B.α) && isapproxinteger(A.β-B.β)
        ct=conversion_type(A.space,B.space)
        ct==NoSpace()?NoSpace():JacobiWeight(max(A.α,B.α),max(A.β,B.β),ct)
    else
        NoSpace()
    end
end

conversion_rule{T,D<:IntervalDomain}(A::JacobiWeight,B::UnivariateSpace{T,D})=conversion_type(A,JacobiWeight(0,0,B))



function Conversion{JS1,JS2,DD<:IntervalDomain}(A::JacobiWeight{JS1,DD},B::JacobiWeight{JS2,DD})
    @assert isapproxinteger(A.α-B.α) && isapproxinteger(A.β-B.β)

    if isapprox(A.α,B.α) && isapprox(A.β,B.β)
        ConversionWrapper(SpaceOperator(Conversion(A.space,B.space),A,B))
    elseif A.space==B.space
        @assert A.α≥B.α&&A.β≥B.β
        d=domain(A)
        x=Fun(identity,d)
        M=tocanonical(d,x)
        m=(1+M).^round(Int,A.α-B.α).*(1-M).^round(Int,A.β-B.β)
        MC=Multiplication(m,B.space)
        # The following is just a safety check
        @assert rangespace(MC) == B.space
        ConversionWrapper(SpaceOperator(MC,A,B))# Wrap the operator with the correct spaces
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
        ConversionWrapper(SpaceOperator(MC,A,B))
    end
end

Conversion{JS,D<:IntervalDomain}(A::RealUnivariateSpace{D},B::JacobiWeight{JS,D})=ConversionWrapper(
    SpaceOperator(
        Conversion(JacobiWeight(0,0,A),B),
        A,B))
Conversion{JS,D<:IntervalDomain}(A::JacobiWeight{JS,D},B::RealUnivariateSpace{D})=ConversionWrapper(
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
                @assert isa(d,IntervalDomain)
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
                @assert isa(d,IntervalDomain)
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

        function getindex{λ,D<:Interval,T}(Σ::$Func{JacobiWeight{Ultraspherical{λ,D},D},T},kr::Range)
            dsp = domainspace(Σ)
            d = domain(Σ)

            @assert dsp.α==dsp.β==λ-0.5

            C = $Len(d)/2

            promote_type(T,typeof(C))[k == 1? C*gamma(λ+one(T)/2)*gamma(one(T)/2)/gamma(λ+one(T)) : zero(T) for k=kr]
        end

        datalength{λ,D<:Interval}(Σ::$Func{JacobiWeight{Ultraspherical{λ,D},D}})=1
    end
end
