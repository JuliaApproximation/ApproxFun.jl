##  Jacobi Operator




recA{T,λ}(::Type{T},::Ultraspherical{λ},k) = (2*(k+λ))/(k+one(T))   # one(T) ensures we get correct type
recB{T}(::Type{T},::Ultraspherical,::) = zero(T)
recC{T,λ}(::Type{T},::Ultraspherical{λ},k) = (k-one(T)+2λ)/(k+one(T))   # one(T) ensures we get correct type

# x p_k
recα{T}(::Type{T},::Ultraspherical,::) = zero(T)
recβ{T,λ}(::Type{T},::Ultraspherical{λ},k) = k/(2*(k-one(T)+λ))   # one(T) ensures we get correct type
recγ{T,λ}(::Type{T},::Ultraspherical{λ},k) = (k-2+2λ)/(2*(k-one(T)+λ))   # one(T) ensures we get correct type



## Multiplication
# these are special cases


Base.stride{U<:Chebyshev,V<:Ultraspherical}(M::ConcreteMultiplication{U,V}) =
    stride(M.f)
Base.stride{U<:Ultraspherical,V<:Chebyshev}(M::ConcreteMultiplication{U,V}) =
    stride(M.f)
Base.stride{U<:Ultraspherical,V<:Ultraspherical}(M::ConcreteMultiplication{U,V}) =
    stride(M.f)




function getindex{D,T,C<:Chebyshev}(M::ConcreteMultiplication{C,Ultraspherical{1,D},T},k::Integer,j::Integer)
    cfs=coefficients(M.f)
    toeplitz_getindex(.5cfs,k,j)+hankel_getindex(-.5cfs[3:end],k,j)
end



function Base.convert{C<:Chebyshev,D,V,T}(::Type{BandedMatrix},S::SubOperator{T,ConcreteMultiplication{C,Ultraspherical{1,D},V,T},Tuple{UnitRange{Int},UnitRange{Int}}})
    ret=bzeros(S)

    kr,jr=parentindexes(S)
    cfs=parent(S).f.coefficients

    # Toeplitz part
    sym_toeplitz_axpy!(1.0,0.5,cfs,kr,jr,ret)

    #Hankel part
    hankel_axpy!(-0.5,@compat(view(cfs,3:length(cfs))),kr,jr,ret)

    ret
end






## Derivative


#Derivative(k::Integer,d::IntervalDomain)=Derivative(k-1:k,d)
#Derivative(d::IntervalDomain)=Derivative(1,d)


Derivative{λ,DD<:Interval}(sp::Ultraspherical{λ,DD},order::Integer) =
    ConcreteDerivative(sp,order)
Integral{λ,DD<:Interval}(sp::Ultraspherical{λ,DD},order::Integer) =
    ConcreteIntegral(sp,order)


rangespace{λ,DD<:Interval}(D::ConcreteDerivative{Ultraspherical{λ,DD}}) =
    Ultraspherical{λ+D.order}(domain(D))
bandinds{λ,DD<:Interval}(D::ConcreteDerivative{Ultraspherical{λ,DD}}) = 0,D.order
bandinds{λ,DD<:Interval}(D::ConcreteIntegral{Ultraspherical{λ,DD}}) = -D.order,0
Base.stride{λ,DD<:Interval}(D::ConcreteDerivative{Ultraspherical{λ,DD}}) = D.order


function getindex{λ,DD<:Interval,K,T}(D::ConcreteDerivative{Ultraspherical{λ,DD},K,T},k::Integer,j::Integer)
    m=D.order
    d=domain(D)

    if j==k+m
        (pochhammer(one(T)*λ,m)*(4./(d.b-d.a)).^m)::T
    else
        zero(T)
    end
end


## Integral

linesum{λ,DD<:Interval}(f::Fun{Ultraspherical{λ,DD}}) =
    sum(setcanonicaldomain(f))*arclength(d)/2





rangespace{λ,DD<:Interval}(D::ConcreteIntegral{Ultraspherical{λ,DD}}) =
    λ == 1 ? Chebyshev() : Ultraspherical{λ-D.order}(domain(D))

function getindex{λ,DD<:Interval,T}(D::ConcreteIntegral{Ultraspherical{λ,DD},T},k::Integer,j::Integer)
    m=D.order
    d=domain(D)
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


function Conversion{b,DD}(A::Chebyshev,B::Ultraspherical{b,DD})
    if b==1
        ConcreteConversion(A,B)
    else
        d=domain(A)
        ConversionWrapper(TimesOperator(Conversion(Ultraspherical{1,DD}(d),B),
                                        Conversion(Chebyshev(d),Ultraspherical{1,DD}(d))))
    end
end


function Conversion{a,b,DD}(A::Ultraspherical{a,DD},B::Ultraspherical{b,DD})
    @assert b >= a

    if b==a
        ConversionWrapper(eye(A))
    elseif b==a+1
        ConcreteConversion(A,B)
    else
        d=domain(A)
        ConversionWrapper(TimesOperator(Conversion(Ultraspherical{b-1,DD}(d),B),
                                        Conversion(A,Ultraspherical{b-1,DD}(d))))
    end
end


function getindex{DD,C<:Chebyshev,T}(M::ConcreteConversion{C,Ultraspherical{1,DD},T},
                                     k::Integer,j::Integer)
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


function getindex{m,λ,DD,T}(M::ConcreteConversion{Ultraspherical{m,DD},Ultraspherical{λ,DD},T},
                            k::Integer,j::Integer)
    #  we can assume that λ==m+1
    c=λ-one(T)  # this supports big types
    if k==j
        c/(k - 2 + λ)
    elseif j==k+2
        -c/(k + λ)
    else
        zero(T)
    end
end


bandinds{DD}(C::ConcreteConversion{Chebyshev{DD},Ultraspherical{1,DD}}) = 0,2
bandinds{m,DD,λ}(C::ConcreteConversion{Ultraspherical{m,DD},Ultraspherical{λ,DD}}) =
    0,2

Base.stride{λ,DD}(C::ConcreteConversion{Chebyshev{DD},Ultraspherical{λ,DD}}) = 2
Base.stride{m,λ,DD}(C::ConcreteConversion{Ultraspherical{m,DD},Ultraspherical{λ,DD}}) = 2


## coefficients

# return the space that has banded Conversion to the other
conversion_rule{border}(a::Chebyshev,b::Ultraspherical{border}) =
    if domainscompatible(a,b)
        a
    else
        NoSpace()
    end

conversion_rule{aorder,border}(a::Ultraspherical{aorder},b::Ultraspherical{border}) =
    if domainscompatible(a,b)
        aorder < border?a:b
    else
        NoSpace()
    end



coefficients(g::Vector,::Ultraspherical{1},::Chebyshev)=ultraiconversion(g)
coefficients(g::Vector,::Chebyshev,::Ultraspherical{1})=ultraconversion(g)


# TODO: include in getindex to speed up
Integral{DD<:Interval}(sp::Chebyshev{DD},m::Integer) =
    IntegralWrapper(TimesOperator([Integral(Ultraspherical{m}(domain(sp)),m),
                                   Conversion(sp,Ultraspherical{m}(domain(sp)))]),m)
