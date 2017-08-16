##  Jacobi Operator




recA(::Type{T},S::Ultraspherical,k) where {T} = (2*(k+order(S)))/(k+one(T))   # one(T) ensures we get correct type
recB(::Type{T},::Ultraspherical,::) where {T} = zero(T)
recC(::Type{T},S::Ultraspherical,k) where {T} = (k-one(T)+2order(S))/(k+one(T))   # one(T) ensures we get correct type

# x p_k
recα(::Type{T},::Ultraspherical,::) where {T} = zero(T)
recβ(::Type{T},S::Ultraspherical,k) where {T} = k/(2*(k-one(T)+order(S)))   # one(T) ensures we get correct type
recγ(::Type{T},S::Ultraspherical,k) where {T} = (k-2+2order(S))/(2*(k-one(T)+order(S)))   # one(T) ensures we get correct type



## Multiplication
# these are special cases


Base.stride(M::ConcreteMultiplication{U,V}) where {U<:Chebyshev,V<:Ultraspherical} =
    stride(M.f)
Base.stride(M::ConcreteMultiplication{U,V}) where {U<:Ultraspherical,V<:Chebyshev} =
    stride(M.f)
Base.stride(M::ConcreteMultiplication{U,V}) where {U<:Ultraspherical,V<:Ultraspherical} =
    stride(M.f)


function Multiplication(f::Fun{C},sp::Ultraspherical{Int}) where C<:Chebyshev
    if order(sp) == 1
        cfs = f.coefficients
        MultiplicationWrapper(f,
            SpaceOperator(SymToeplitzOperator(cfs/2) +
                                HankelOperator(view(cfs,3:length(cfs))/(-2)),
                          sp,sp))

    else
        ConcreteMultiplication(f,sp)
    end
end


## Derivative


#Derivative(k::Integer,d::IntervalDomain)=Derivative(k-1:k,d)
#Derivative(d::IntervalDomain)=Derivative(1,d)


Derivative(sp::Ultraspherical{LT,DD},m::Integer) where {LT,DD<:Segment} =
    ConcreteDerivative(sp,m)
function Integral(sp::Ultraspherical{LT,DD},m::Integer) where {LT,DD<:Segment}
    λ = order(sp)
    if m ≤ λ
        ConcreteIntegral(sp,m)
    else # Convert up
        nsp = Ultraspherical(λ+1,domain(sp))
        IntegralWrapper(Integral(nsp,m)*Conversion(sp,nsp),m)
    end
end


rangespace(D::ConcreteDerivative{Ultraspherical{LT,DD,RR}}) where {LT,DD<:Segment,RR} =
    Ultraspherical(order(domainspace(D))+D.order,domain(D))
bandinds(D::ConcreteDerivative{Ultraspherical{LT,DD,RR}}) where {LT,DD<:Segment,RR} = D.order,D.order
bandinds(D::ConcreteIntegral{Ultraspherical{LT,DD,RR}}) where {LT,DD<:Segment,RR} = -D.order,-D.order
Base.stride(D::ConcreteDerivative{Ultraspherical{LT,DD,RR}}) where {LT,DD<:Segment,RR} = D.order


function getindex(D::ConcreteDerivative{Ultraspherical{TT,DD,RR},K,T},
               k::Integer,j::Integer) where {TT,DD<:Segment,RR,K,T}
    m=D.order
    d=domain(D)
    λ=order(domainspace(D))

    if j==k+m
        T((pochhammer(one(T)*λ,m)*(4./(d.b-d.a)).^m))
    else
        zero(T)
    end
end


## Integral

linesum(f::Fun{Ultraspherical{LT,DD,RR}}) where {LT,DD<:Segment,RR} =
    sum(setcanonicaldomain(f))*arclength(d)/2





rangespace(D::ConcreteIntegral{Ultraspherical{LT,DD,RR}}) where {LT,DD<:Segment,RR} =
    order(domainspace(D)) == 1 ? Chebyshev() : Ultraspherical(order(domainspace(D))-D.order,domain(D))

function getindex(Q::ConcreteIntegral{Ultraspherical{LT,DD,RR}},k::Integer,j::Integer) where {LT,DD<:Segment,RR}
    T=eltype(Q)
    m=Q.order
    d=domain(Q)
    λ=order(domainspace(Q))
    @assert m<=λ

    if λ == 1 && k==j+1
        C = (d.b-d.a)/2
        T(C./(k-1))
    elseif λ > 1 && k==j+m
        T(pochhammer(one(T)*λ,-m)*((d.b-d.a)/4)^m)
    else
        zero(T)
    end
end



## Conversion Operator


function Conversion(A::Chebyshev,B::Ultraspherical)
    if order(B) ≤ 1
        ConcreteConversion(A,B)
    else
        d=domain(A)
        US=Ultraspherical(order(B)-1,d)
        ConversionWrapper(TimesOperator(Conversion(US,B),
                                        Conversion(Chebyshev(d),US)))
    end
end

function Conversion(A::Ultraspherical,B::Chebyshev)
    if order(A) == 1//2
        ConcreteConversion(A,B)
    else
        error("Not implemented")
    end
end


maxspace_rule(A::Ultraspherical,B::Chebyshev) = A


function Conversion(A::Ultraspherical,B::Ultraspherical)
    a=order(A); b=order(B)
    if b==a
        ConversionWrapper(eye(A))
    elseif a<b≤a+1  || b<a≤b+1
        ConcreteConversion(A,B)
    elseif b ≠ 1
        d=domain(A)
        ConversionWrapper(TimesOperator(Conversion(Ultraspherical(b-1,d),B),
                                        Conversion(A,Ultraspherical(b-1,d))))
    else
        error("Cannot convert from $A to $B")
    end
end

maxspace_rule(A::Ultraspherical,B::Ultraspherical) = order(A) > order(B) ? A : B


function getindex(M::ConcreteConversion{C,Ultraspherical{Int,DD,RR},T},
               k::Integer,j::Integer) where {DD,RR,C<:Chebyshev,T}
   # order must be 1
    if k==j==1
        one(T)
    elseif k==j
        one(T)/2
    elseif j==k+2
        -one(T)/2
    else
        zero(T)
    end
end


function getindex(M::ConcreteConversion{Ultraspherical{Int,DD,RR},Ultraspherical{Int,DD,RR},T},
                   k::Integer,j::Integer) where {DD,RR,T}
    #  we can assume that λ==m+1
    λ=order(rangespace(M))
    c=λ-one(T)  # this supports big types
    if k==j
        c/(k - 2 + λ)
    elseif j==k+2
        -c/(k + λ)
    else
        zero(T)
    end
end

function getindex(M::ConcreteConversion{Ultraspherical{LT,DD,RR},Ultraspherical{LT,DD,RR},T},
                k::Integer,j::Integer) where {LT,DD,RR,T}
    λ=order(rangespace(M))
    if order(domainspace(M))+1==λ
        c=λ-one(T)  # this supports big types
        if k==j
            c/(k - 2 + λ)
        elseif j==k+2
            -c/(k + λ)
        else
            zero(T)
        end
    else
        error("Not implemented")
    end
end


bandinds(C::ConcreteConversion{Chebyshev{DD,RR},Ultraspherical{Int,DD,RR}}) where {DD,RR} = 0,2  # order == 1
bandinds(C::ConcreteConversion{Ultraspherical{Int,DD,RR},Ultraspherical{Int,DD,RR}}) where {DD,RR} = 0,2

bandinds(C::ConcreteConversion{Chebyshev{DD,RR},Ultraspherical{LT,DD,RR}}) where {LT,DD,RR} =
    0,order(rangespace(C))==1?2:∞
bandinds(C::ConcreteConversion{Ultraspherical{LT,DD,RR},Chebyshev{DD,RR}}) where {LT,DD,RR} =
    0,order(domainspace(C))==1?2:∞

bandinds(C::ConcreteConversion{Ultraspherical{LT1,DD,RR},Ultraspherical{LT2,DD,RR}}) where {LT1,LT2,DD,RR} =
    0,order(domainspace(C))+1==order(rangespace(C))?2:∞

Base.stride(C::ConcreteConversion{Chebyshev{DD,RR},Ultraspherical{Int,DD,RR}}) where {DD,RR} = 2
Base.stride(C::ConcreteConversion{Ultraspherical{LT1,DD,RR},Ultraspherical{LT2,DD,RR}}) where {LT1,LT2,DD,RR} = 2


## coefficients

# return the space that has banded Conversion to the other
conversion_rule(a::Chebyshev,b::Ultraspherical{Int}) =
    if domainscompatible(a,b)
        a
    else
        NoSpace()
    end

conversion_rule(a::Ultraspherical{LT},b::Ultraspherical{LT}) where {LT} =
    if domainscompatible(a,b) && isapproxinteger(order(a)-order(b))
        order(a) < order(b)?a:b
    else
        NoSpace()
    end



function coefficients(g::AbstractVector,sp::Ultraspherical{Int},C::Chebyshev)
    if order(sp) == 1
        ultraiconversion(g)
    else
        # do one at a time
        coefficients(g,sp,Ultraspherical(1,domain(sp)),C)
    end
end
function coefficients(g::AbstractVector,C::Chebyshev,sp::Ultraspherical)
    if order(sp) == 1
        ultraconversion(g)
    else
        # do one at a time
        coefficients(g,C,Ultraspherical(1,domain(sp)),sp)
    end
end


# TODO: include in getindex to speed up
Integral(sp::Chebyshev{DD},m::Integer) where {DD<:Segment} =
    IntegralWrapper(TimesOperator([Integral(Ultraspherical(m,domain(sp)),m),
                                   Conversion(sp,Ultraspherical(m,domain(sp)))]),m)




## Non-banded conversions

function getindex(M::ConcreteConversion{C,Ultraspherical{LT,DD,RR},T},
            k::Integer,j::Integer) where {DD,RR,LT,C<:Chebyshev,T}
    λ = order(rangespace(M))
    if λ == 1
        if k==j==1
            one(T)
        elseif k==j
            one(T)/2
        elseif j==k+2
            -one(T)/2
        else
            zero(T)
        end
    elseif λ == 0.5
        # Cheb-to-Leg
        if j==k==1
            one(T)
        elseif j==k
            T(sqrt(π)/(2FastTransforms.Λ(k-1)))
        elseif k < j && iseven(k-j)
            T(-(j-1)*(k-0.5)*(FastTransforms.Λ((j-k-2)/2)/(j-k))*
                            (FastTransforms.Λ((j+k-3)/2)/(j+k-1)))
        else
            zero(T)
        end
    else
        error("Not implemented")
    end
end


function getindex(M::ConcreteConversion{Ultraspherical{LT,DD,RR},C,T},
            k::Integer,j::Integer) where {DD,RR,LT,C<:Chebyshev,T}
    λ = order(domainspace(M))
    if λ == 1
        # order must be 1
        if k==j==1
            one(T)
        elseif k==j
            one(T)/2
        elseif j==k+2
            -one(T)/2
        else
            zero(T)
        end
    elseif λ == 0.5
        if k==1 && isodd(j)
            T(FastTransforms.Λ((j-1)/2)^2/π)
        elseif k ≤ j && iseven(k-j)
            T(FastTransforms.Λ((j-k)/2)*FastTransforms.Λ((k+j-2)/2)*2/π)
        else
            zero(T)
        end
    else
        error("Not implemented")
    end
end



function getindex(M::ConcreteConversion{Ultraspherical{LT,DD,RR},
                                     Ultraspherical{LT2,DD,RR},T},
                     k::Integer,j::Integer) where {DD,RR,LT,LT2,T}
    λ1 = order(domainspace(M))
    λ2 = order(rangespace(M))
    if abs(λ1-λ2) < 1
        if j ≥ k && iseven(k-j)
            T((λ1 < λ2 && k ≠ j ? -1 : 1) *  # fix sign for lgamma
                exp(lgamma(λ2)+log(k-1+λ2)-lgamma(λ1)-lgamma(λ1-λ2) + lgamma((j-k)/2+λ1-λ2)-
                lgamma((j-k)/2+1)+lgamma((k+j-2)/2+λ1)-lgamma((k+j-2)/2+λ2+1)))
        else
            zero(T)
        end
    else
        error("Not implemented")
    end
end


ReverseOrientation(S::Ultraspherical) = ReverseOrientationWrapper(SpaceOperator(NegateEven(),S,reverseorientation(S)))
Reverse(S::Ultraspherical) = ReverseWrapper(SpaceOperator(NegateEven(),S,S))
