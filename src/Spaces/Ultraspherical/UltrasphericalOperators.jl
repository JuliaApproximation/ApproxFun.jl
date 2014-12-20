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
function Base.getindex(op::Evaluation{Chebyshev,Bool},k::Range)
    x = op.x
    d = domain(op)
    p = op.order
    cst = (2/(d.b-d.a))^p

    if x
        ret = ones(size(k)[1])
    elseif !x
        ret = (-1.)^(p+1)*(-1.).^k ##TODO: speed up
    end

    for m=0:p-1
        ret .*= ((k-1).^2-m^2)./(2m+1)
    end

    return ret*cst
end

function Base.getindex(op::Evaluation{Chebyshev},k::Range)
    if op.order == 0    
        evaluatechebyshev(k[end],tocanonical(domain(op),op.x))[k]
    else
        error("Only zero–second order implemented")
    end
end



## Multiplication

function chebmult_addentries!(cfs::Vector,A::ShiftArray,kr::Range)
    toeplitz_addentries!(.5cfs,A,kr)
    hankel_addentries!(.5cfs,A,max(kr[1],2):kr[end])            
end


addentries!{D<:Ultraspherical}(M::Multiplication{D,Chebyshev},A::ShiftArray,kr::Range)=chebmult_addentries!(canonicalcoefficients(M.f),A,kr)

function addentries!{D<:Ultraspherical}(M::Multiplication{D,Ultraspherical{1}},A::ShiftArray,kr::Range)
    cfs=canonicalcoefficients(M.f)
    toeplitz_addentries!(.5cfs,A,kr)
    hankel_addentries!(-.5cfs[3:end],A,kr)    
end



function usjacobi_addentries!(λ::Integer,A::ShiftArray,kr::Range1)
    for k=kr
        A[k,-1]=.5(k-1)/(k-2+λ)
        A[k,1]=.5(k+2λ-1)/(k+λ)
    end
    A
end

function addentries!{D<:Ultraspherical,λ}(M::Multiplication{D,Ultraspherical{λ}},A::ShiftArray,kr::Range)
    a=coefficients(M.f,domainspace(M))
    for k=kr
        A[k,0]=a[1] 
    end

    if length(a) > 1
        jkr=max(1,kr[1]-length(a)+1):kr[end]+length(a)-1
        ##TODO: simplify shift array and combine with Ultraspherical        
        J=BandedArray(ShiftArray(zeros(length(jkr),3),1-jkr[1],2),jkr)
        usjacobi_addentries!(λ,J.data,jkr)
    
        C1=2λ*J
    
        shiftarray_const_addentries!(C1.data,a[2],A,kr)

        C0=BandedArray(ShiftArray(ones(length(jkr),1),1-jkr[1],1),jkr)
    
        for k=1:length(a)-2    
            C1,C0=2(k+λ)./(k+1)*J*C1-(k+2λ-1)./(k+1).*C0,C1
            shiftarray_const_addentries!(C1.data,a[k+2],A,kr)    
        end
    end
    
    A
end



## Derivative


#Derivative(k::Integer,d::IntervalDomain)=Derivative(k-1:k,d)
#Derivative(d::IntervalDomain)=Derivative(1,d)


rangespace{λ}(D::Derivative{Ultraspherical{λ}})=Ultraspherical{λ+D.order}(domain(D))
bandinds{S<:Ultraspherical}(D::Derivative{S})=0,D.order
bandinds{S<:Ultraspherical}(D::Integral{S})=-D.order,0   

function addentries!{λ}(D::Derivative{Ultraspherical{λ}},A::ShiftArray,kr::Range1)
    m=D.order
    d=domain(D)

    @assert isa(d,Interval)

    if λ == 0
        C=2.^(m-1).*factorial(m-1)*(2./(d.b-d.a)).^m    
        for k=kr
            A[k,m] += C*(m+k-1)
        end
    else
        C=2.^m.*factorial(λ+m-1)./factorial(λ-1)*(2./(d.b-d.a)).^m        
        for k=kr        
            A[k,m] += C
        end
    end
    
    A
end
 
## TODO: reimplement

# function Derivative(order::Range1,d::IntervalDomain)
#     @assert order[1] == 0 && order[end] <= 2  ##TODO other orders
#     
#     Mp = Fun(x->tocanonicalD(d,x),d)
#     
#     if order[end] == 1
#         Mp*USDerivative(0:1)
#     elseif order[end] == 2
#         (Mp.^2)*USDerivative(0:2) + diff(Mp)*USDerivative(0:1)
#     end
# end

 
## Integral


rangespace{λ}(D::Integral{Ultraspherical{λ}})=Ultraspherical{λ-D.order}(domain(D))

function addentries!{λ}(D::Integral{Ultraspherical{λ}},A::ShiftArray,kr::Range1)
    m=D.order
    d=domain(D)
    @assert m<=λ
    @assert isa(d,Interval)

    if λ == 1
        C = .5(d.b-d.a)
        for k=max(kr[1],2):kr[end]
            A[k,-1] += C./(k-1)
        end
    else
        C=factorial(λ-m-1)/factorial(λ-1)/2.^m*(.5(d.b-d.a))^m        
        for k=kr
            A[k,-m] += C
        end
    end
    
    A
end



## Conversion Operator




function Conversion{a,b}(A::Ultraspherical{a},B::Ultraspherical{b})
    @assert b > a

    if b==a+1
        Conversion{Ultraspherical{a},Ultraspherical{b},Float64}(A,B)
    else
        d=domain(A)
        Conversion(Ultraspherical{b-1}(d),B)*Conversion(A,Ultraspherical{b-1}(d))
    end
end   


function addentries!(M::Conversion{Chebyshev,Ultraspherical{1}},A::ShiftArray,kr::Range)
    for k=kr
        A[k,0] += (k == 1)? 1. : .5
        A[k,2] += -.5        
    end
    
    A    
end

function addentries!{m,λ}(M::Conversion{Ultraspherical{m},Ultraspherical{λ}},A::ShiftArray,kr::Range)
    @assert λ==m+1
    for k=kr
        A[k,0] += (λ-1.)./(k - 2. + λ)
        A[k,2] += -(λ-1.)./(k + λ)
    end
    
    A    
end

function multiplyentries!(M::Conversion{Chebyshev,Ultraspherical{1}},A::ShiftArray,kr::Range)
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

function multiplyentries!{m,λ}(M::Conversion{Ultraspherical{m},Ultraspherical{λ}},A::ShiftArray,kr::Range)
    @assert λ==m+1
    cr=columnrange(A)::Range1{Int64}
    
    λf = 1.λ
    
    #We assume here that the extra rows are redundant
    for k=max(kr[1],1):kr[end]+2,j=cr
        A[k,j] *= (λf-1)./(k - 2. + λf)
    end
    
    #We assume that A has allocated 2 more bandwidth
    for k=max(kr[1],1):kr[end],j=(cr[1]+2):cr[end]
        A[k,j] -= A[k+2,j-2]
    end 
end

bandinds{m,λ}(C::Conversion{Ultraspherical{m},Ultraspherical{λ}})=0,2



## spaceconversion

# return the space that has banded Conversion to the other
function conversion_rule{aorder,border}(a::Ultraspherical{aorder},b::Ultraspherical{border})
    @assert domainscompatible(a,b)
    
    aorder < border?a:b
end


spaceconversion(g::Vector,::Ultraspherical{1},::Chebyshev)=ultraiconversion(g)
spaceconversion(g::Vector,::Chebyshev,::Ultraspherical{1})=ultraconversion(g)



