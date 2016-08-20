##  Jacobi Operator




recA{T}(::Type{T},S::Ultraspherical,k) = (2*(k+order(S)))/(k+one(T))   # one(T) ensures we get correct type
recB{T}(::Type{T},::Ultraspherical,::) = zero(T)
recC{T}(::Type{T},S::Ultraspherical,k) = (k-one(T)+2order(S))/(k+one(T))   # one(T) ensures we get correct type

# x p_k
recα{T}(::Type{T},::Ultraspherical,::) = zero(T)
recβ{T}(::Type{T},S::Ultraspherical,k) = k/(2*(k-one(T)+order(S)))   # one(T) ensures we get correct type
recγ{T}(::Type{T},S::Ultraspherical,k) = (k-2+2order(S))/(2*(k-one(T)+order(S)))   # one(T) ensures we get correct type



## Multiplication
# these are special cases


Base.stride{U<:Chebyshev,V<:Ultraspherical}(M::ConcreteMultiplication{U,V}) =
    stride(M.f)
Base.stride{U<:Ultraspherical,V<:Chebyshev}(M::ConcreteMultiplication{U,V}) =
    stride(M.f)
Base.stride{U<:Ultraspherical,V<:Ultraspherical}(M::ConcreteMultiplication{U,V}) =
    stride(M.f)


function Multiplication{C<:Chebyshev}(f::Fun{C},sp::Ultraspherical{Int})
    if order(sp) == 1
        cfs = f.coefficients
        MultiplicationWrapper(f,
            SpaceOperator(SymToeplitzOperator(cfs/2) +
                                HankelOperator(@compat(view(cfs,3:length(cfs)))/(-2)),
                          sp,sp))

    else
        ConcreteMultiplication(f,sp)
    end
end


## Derivative


#Derivative(k::Integer,d::IntervalDomain)=Derivative(k-1:k,d)
#Derivative(d::IntervalDomain)=Derivative(1,d)


Derivative{LT,DD<:Interval}(sp::Ultraspherical{LT,DD},order::Integer) =
    ConcreteDerivative(sp,order)
Integral{LT,DD<:Interval}(sp::Ultraspherical{LT,DD},order::Integer) =
    ConcreteIntegral(sp,order)


rangespace{LT,DD<:Interval}(D::ConcreteDerivative{Ultraspherical{LT,DD}}) =
    Ultraspherical(order(domainspace(D))+D.order,domain(D))
bandinds{LT,DD<:Interval}(D::ConcreteDerivative{Ultraspherical{LT,DD}}) = 0,D.order
bandinds{LT,DD<:Interval}(D::ConcreteIntegral{Ultraspherical{LT,DD}}) = -D.order,0
Base.stride{LT,DD<:Interval}(D::ConcreteDerivative{Ultraspherical{LT,DD}}) = D.order


function getindex{TT,DD<:Interval,K,T}(D::ConcreteDerivative{Ultraspherical{TT,DD},K,T},
                                      k::Integer,j::Integer)
    m=D.order
    d=domain(D)
    λ=order(domainspace(D))

    if j==k+m
        (pochhammer(one(T)*λ,m)*(4./(d.b-d.a)).^m)::T
    else
        zero(T)
    end
end


## Integral

linesum{LT,DD<:Interval}(f::Fun{Ultraspherical{LT,DD}}) =
    sum(setcanonicaldomain(f))*arclength(d)/2





rangespace{LT,DD<:Interval}(D::ConcreteIntegral{Ultraspherical{LT,DD}}) =
    order(domainspace(D)) == 1 ? Chebyshev() : Ultraspherical(order(domainspace(D))-D.order,domain(D))

function getindex{LT,DD<:Interval,T}(D::ConcreteIntegral{Ultraspherical{LT,DD},T},k::Integer,j::Integer)
    m=D.order
    d=domain(D)
    λ=order(domainspace(D))
    @assert m<=λ

    if λ == 1 && k==j+1
        C = .5(d.b-d.a)
        C./(k-one(T))
    elseif λ > 1 && k==j+m
        pochhammer(one(T)*λ,-m)*(.25(d.b-d.a))^m
    else
        zero(T)
    end
end



## Conversion Operator


function Conversion{DD}(A::Chebyshev,B::Ultraspherical{Int,DD})
    if order(B) == 1
        ConcreteConversion(A,B)
    else
        d=domain(A)
        ConversionWrapper(TimesOperator(Conversion(Ultraspherical(1,d),B),
                                        Conversion(Chebyshev(d),Ultraspherical(1,d))))
    end
end


function Conversion{LT,DD}(A::Ultraspherical{LT,DD},B::Ultraspherical{LT,DD})
    a=order(A); b=order(B)
    @assert b >= a

    if b==a
        ConversionWrapper(eye(A))
    elseif b==a+1
        ConcreteConversion(A,B)
    else
        d=domain(A)
        ConversionWrapper(TimesOperator(Conversion(Ultraspherical(b-1,d),B),
                                        Conversion(A,Ultraspherical(b-1,d))))
    end
end


function getindex{DD,C<:Chebyshev,T}(M::ConcreteConversion{C,Ultraspherical{Int,DD},T},
                                     k::Integer,j::Integer)
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


function getindex{LT,DD,T}(M::ConcreteConversion{Ultraspherical{LT,DD},Ultraspherical{LT,DD},T},
                            k::Integer,j::Integer)
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


bandinds{DD}(C::ConcreteConversion{Chebyshev{DD},Ultraspherical{Int,DD}}) = 0,2  # order == 1
bandinds{LT,DD}(C::ConcreteConversion{Ultraspherical{LT,DD},Ultraspherical{LT,DD}}) =
    0,2

Base.stride{DD}(C::ConcreteConversion{Chebyshev{DD},Ultraspherical{Int,DD}}) = 2
Base.stride{LT,DD}(C::ConcreteConversion{Ultraspherical{LT,DD},Ultraspherical{LT,DD}}) = 2


## coefficients

# return the space that has banded Conversion to the other
conversion_rule(a::Chebyshev,b::Ultraspherical{Int}) =
    if domainscompatible(a,b)
        a
    else
        NoSpace()
    end

conversion_rule{LT}(a::Ultraspherical{LT},b::Ultraspherical{LT}) =
    if domainscompatible(a,b) && isapproxinteger(order(a)-order(b))
        order(a) < order(b)?a:b
    else
        NoSpace()
    end



function coefficients(g::Vector,sp::Ultraspherical{Int},C::Chebyshev)
    if order(sp) == 1
        ultraiconversion(g)
    else
        # do one at a time
        coefficients(g,sp,Ultraspherical(1,domain(sp)),C)
    end
end
function coefficients(g::Vector,C::Chebyshev,sp::Ultraspherical)
    if order(sp) == 1
        ultraconversion(g)
    else
        # do one at a time
        coefficients(g,C,Ultraspherical(1,domain(sp)),sp)
    end
end


# TODO: include in getindex to speed up
Integral{DD<:Interval}(sp::Chebyshev{DD},m::Integer) =
    IntegralWrapper(TimesOperator([Integral(Ultraspherical(m,domain(sp)),m),
                                   Conversion(sp,Ultraspherical(m,domain(sp)))]),m)
