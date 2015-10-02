##  Jacobi Operator

recB{T}(::Type{T},::Ultraspherical,::)=zero(T)

recA{T}(::Type{T},::Chebyshev,k)=2one(T)
recC{T}(::Type{T},::Chebyshev,k)=one(T)   # one(T) ensures we get correct type

recA{T,λ}(::Type{T},::Ultraspherical{λ},k)=(2*(k+λ))/(k+one(T))   # one(T) ensures we get correct type
recC{T,λ}(::Type{T},::Ultraspherical{λ},k)=(k-one(T)+2λ)/(k+one(T))   # one(T) ensures we get correct type

# x p_k
recα{T}(::Type{T},::Ultraspherical,::)=zero(T)

recβ{T}(::Type{T},::Chebyshev,k)=ifelse(k==1,one(T),one(T)/2)   # one(T) ensures we get correct type,ifelse ensures inlining
recγ{T}(::Type{T},::Chebyshev,k)=one(T)/2   # one(T) ensures we get correct type

recβ{T,λ}(::Type{T},::Ultraspherical{λ},k)=k/(2*(k-one(T)+λ))   # one(T) ensures we get correct type
recγ{T,λ}(::Type{T},::Ultraspherical{λ},k)=(k-2+2λ)/(2*(k-one(T)+λ))   # one(T) ensures we get correct type



## Evaluation

function evaluatechebyshev{T<:Number}(n::Integer,x::T)
    if n == 1
        [one(T)]
    else
        p = zeros(T,n)
        p[1] = one(T)
        p[2] = x

        for j=2:n-1
            p[j+1] = 2x*p[j] - p[j-1]
        end

        p
    end
end





##TODO: the overloading as both vector and row vector may be confusing
function Base.getindex{T<:Number,C<:Chebyshev}(op::Evaluation{C,Bool,T},k::Range)
    x = op.x
    d = domain(op)
    p = op.order
    cst = (2/(d.b-d.a))^p

    if x
        ret = ones(T,length(k))
    else
        ret = Array(T,length(k))
        k1=1-first(k)
        for j=k
            ret[j+k1]=(-1)^(p+1)*(-one(T))^j
        end
    end

    for m=0:p-1
        k1=1-first(k)
        for j=k
            ret[j+k1] *= (j-1)^2-m^2
        end
        ret /= 2m+1
    end

    return ret*cst
end

function Base.getindex{C<:Chebyshev}(op::Evaluation{C},k::Range)
    if op.order == 0
        evaluatechebyshev(k[end],tocanonical(domain(op),op.x))[k]
    else
        error("Only zero–second order implemented")
    end
end



## Multiplication
# these are special cases


Base.stride{U<:Ultraspherical,V<:Ultraspherical}(M::Multiplication{U,V})=stride(M.f)


function chebmult_addentries!(cfs::Vector,A,kr::Range)
    toeplitz_addentries!(.5cfs,A,kr)
    hankel_addentries!(.5cfs,A,intersect(2:kr[end],kr))
end


for TYP in (:UnitRange,:Range) # needed to avoid confusion
    @eval begin
        addentries!{T,C<:Chebyshev}(M::Multiplication{C,C,T},A,kr::$TYP,::Colon)=chebmult_addentries!(coefficients(M.f),A,kr)

        function addentries!{D,T,C<:Chebyshev}(M::Multiplication{C,Ultraspherical{1,D},T},A,kr::$TYP,::Colon)
            cfs=coefficients(M.f)
            toeplitz_addentries!(.5cfs,A,kr)
            hankel_addentries!(-.5cfs[3:end],A,kr)
        end
    end
end

function addentries!{λ,PS<:PolynomialSpace,D,T}(M::Multiplication{Ultraspherical{λ,D},PS,T},A,kr::UnitRange,::Colon)
    a=coefficients(M.f)
    for k=kr
        A[k,k]=a[1]
    end

    if length(a) > 1
        jkr=max(1,kr[1]-length(a)+1):kr[end]+length(a)-1

        J=subview(Recurrence(domainspace(M)),jkr,jkr)
        C1=2λ*J
        addentries!(C1,a[2],A,kr,:)
        C0=isbaeye(jkr)

        for k=1:length(a)-2
            C1,C0=2(k+λ)/(k+one(T))*J*C1-(k+2λ-one(T))/(k+one(T))*C0,C1
            addentries!(C1,a[k+2],A,kr,:)
        end
    end

    A
end


function addentries!{PS<:PolynomialSpace,T,C<:Chebyshev}(M::Multiplication{C,PS,T},A,kr::UnitRange,::Colon)
    a=coefficients(M.f)

    for k=kr
        A[k,k]=a[1]
    end

    if length(M.f) > 1
        sp=M.space
        jkr=max(1,kr[1]-length(a)+1):kr[end]+length(a)-1

        #Multiplication is transpose
        J=subview(Recurrence(sp),jkr,jkr)
        C1=J
        addentries!(C1,a[2],A,kr,:)
        C0=isbaeye(jkr)

        for k=1:length(a)-2
            C1,C0=2J*C1-C0,C1
            addentries!(C1,a[k+2],A,kr,:)
        end
    end

    A
end





## Derivative


#Derivative(k::Integer,d::IntervalDomain)=Derivative(k-1:k,d)
#Derivative(d::IntervalDomain)=Derivative(1,d)


rangespace{λ,DD<:Interval}(D::Derivative{Ultraspherical{λ,DD}})=Ultraspherical{λ+D.order}(domain(D))
bandinds{λ,DD<:Interval}(D::Derivative{Ultraspherical{λ,DD}})=0,D.order
bandinds{λ,DD<:Interval}(D::Integral{Ultraspherical{λ,DD}})=-D.order,0
Base.stride{λ,DD<:Interval}(D::Derivative{Ultraspherical{λ,DD}})=D.order

function addentries!{λ,DD<:Interval}(D::Derivative{Ultraspherical{λ,DD}},A,kr::Range,::Colon)
    m=D.order
    d=domain(D)

    if λ == 0
        C=.5pochhammer(1.,m-1)*(4./(d.b-d.a)).^m
        for k=kr
            A[k,k+m] += C*(m+k-1)
        end
    else
        C=pochhammer(1.λ,m)*(4./(d.b-d.a)).^m
        for k=kr
            A[k,k+m] += C
        end
    end

    A
end


## Integral

function linesum{λ,DD<:Interval}(f::Fun{Ultraspherical{λ,DD}})
    d=domain(f)
    sum(Fun(f.coefficients,Ultraspherical{λ}()))*length(d)/2
end



# TODO: include in addentries! to speed up
Integral{DD<:Interval}(sp::Chebyshev{DD},m::Integer)=IntegralWrapper(
    TimesOperator([Integral(Ultraspherical{m}(domain(sp)),m),Conversion(sp,Ultraspherical{m}(domain(sp)))]),m)

rangespace{λ,DD<:Interval}(D::Integral{Ultraspherical{λ,DD}})=Ultraspherical{λ-D.order}(domain(D))

function addentries!{λ,DD<:Interval}(D::Integral{Ultraspherical{λ,DD}},A,kr::Range,::Colon)
    m=D.order
    d=domain(D)
    @assert m<=λ

    if λ == 1
        C = .5(d.b-d.a)
        for k=max(kr[1],2):kr[end]
            A[k,k-1] += C./(k-1)
        end
    elseif λ > 1
        C=pochhammer(1.λ,-m)*(.25(d.b-d.a))^m
        for k=kr
            A[k,k-m] += C
        end
    end

    A
end



## Conversion Operator




function Conversion{a,b,DD}(A::Ultraspherical{a,DD},B::Ultraspherical{b,DD})
    @assert b >= a

    if b==a
        eye(A)
    elseif b==a+1
        Conversion{Ultraspherical{a,DD},Ultraspherical{b,DD},promote_type(Float64,real(eltype(domain(A))),real(eltype(domain(B))))}(A,B)
    else
        d=domain(A)
        TimesOperator(Conversion(Ultraspherical{b-1,DD}(d),B),Conversion(A,Ultraspherical{b-1,DD}(d)))
    end
end


function addentries!{DD,C<:Chebyshev}(M::Conversion{C,Ultraspherical{1,DD}},A,kr::Range,::Colon)
    # this uses that 0.5 is exact, so no need for special bigfloat def
    for k=kr
        A[k,k] += (k == 1)? 1. : .5
        A[k,k+2] += -.5
    end

    A
end

function addentries!{m,λ,DD,T}(M::Conversion{Ultraspherical{m,DD},Ultraspherical{λ,DD},T},A,kr::Range,::Colon)
    @assert λ==m+1
    c=λ-one(T)  # this supports big types
    for k=kr
        A[k,k] += c/(k - 2 + λ)
        A[k,k+2] += -c/(k + λ)
    end

    A
end

function multiplyentries!{DD,C<:Chebyshev}(M::Conversion{C,Ultraspherical{1,DD}},A,kr::Range)
    cr=columnrange(A)::Range1{Int}

    #We assume here that the extra rows are redundant
    for k=max(2,kr[1]):kr[end]+2,j=cr
        A[k,j] *= .5
    end

    #We assume that A has allocated 2 more bandwidth
    for k=max(1,kr[1]):kr[end],j=(cr[1]+2):cr[end]
        A[k,j] -= A[k+2,j-2]
    end
end

function multiplyentries!{m,λ,DD,T}(M::Conversion{Ultraspherical{m,DD},Ultraspherical{λ,DD},T},A,kr::Range)
    @assert λ==m+1
    cr=columnrange(A)::Range1{Int64}

    c = λ-one(T)
    #We assume here that the extra rows are redundant
    for k=max(kr[1],1):kr[end]+2,j=cr
        A[k,j] *= c./(k - 2. + λ)
    end

    #We assume that A has allocated 2 more bandwidth
    for k=max(kr[1],1):kr[end],j=(cr[1]+2):cr[end]
        A[k,j] -= A[k+2,j-2]
    end
end

bandinds{m,DD,λ}(C::Conversion{Ultraspherical{m,DD},Ultraspherical{λ,DD}})=0,2
Base.stride{m,λ,DD}(C::Conversion{Ultraspherical{m,DD},Ultraspherical{λ,DD}})=2


## coefficients

# return the space that has banded Conversion to the other
function conversion_rule{aorder,border}(a::Ultraspherical{aorder},b::Ultraspherical{border})
    if domainscompatible(a,b)
        aorder < border?a:b
    else
        NoSpace()
    end
end


coefficients(g::Vector,::Ultraspherical{1},::Chebyshev)=ultraiconversion(g)
coefficients(g::Vector,::Chebyshev,::Ultraspherical{1})=ultraconversion(g)
