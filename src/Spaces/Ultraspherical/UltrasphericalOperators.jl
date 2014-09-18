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
function Base.getindex(op::Evaluation{ChebyshevSpace,Bool},k::Range)
    x = op.x
    d = domain(op)
    
    if     !x && op.order ==0
        -(-1.).^k  ##TODO: speed up
    elseif  x && op.order ==0
        ones(size(k)[1])
    elseif !x && op.order ==1
        (k-1).*(k-1).*(-1.).^k*2/(d.b-d.a) 
    elseif  x && op.order ==1
        (k-1).*(k-1)*2/(d.b-d.a) 
    elseif !x && op.order ==2
        -(k.-1).^2.*((k.-1).^2.-1)/3.*(-1.).^k*(2/(d.b-d.a))^2
    elseif  x && op.order ==2
        (k.-1).^2.*((k.-1).^2.-1)/3*(2/(d.b-d.a))^2  
    else
        error("Only zero–second order implemented")
    end
end

function Base.getindex(op::Evaluation{ChebyshevSpace,Float64},k::Range)
    if op.order == 0    
        evaluatechebyshev(k[end],tocanonical(domain(op),op.x))[k]
    else
        error("Only zero–second order implemented")
    end
end



## Multiplication


function addentries!{T,D}(M::Multiplication{T,D,ChebyshevSpace},A::ShiftArray,kr::Range1)
    cfs=canonicalcoefficients(M.f)
    toeplitz_addentries!(.5cfs,A,kr)
    hankel_addentries!(.5cfs,A,max(kr[1],2):kr[end])            
end

function addentries!{T,D}(M::Multiplication{T,D,UltrasphericalSpace{1}},A::ShiftArray,kr::Range1)
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

function addentries!{T,D,λ}(M::Multiplication{T,D,UltrasphericalSpace{λ}},A::ShiftArray,kr::Range)
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


rangespace{λ}(D::Derivative{UltrasphericalSpace{λ}})=UltrasphericalSpace{λ+D.order}(domain(D))

function addentries!{λ}(D::Derivative{UltrasphericalSpace{λ}},A::ShiftArray,kr::Range1)
    m=D.order
    d=domain(D)
    μ=λ+m

    @assert isa(d,Interval)

    if λ == 0
        C=2.^(m-1).*factorial(μ-1)*(2./length(d)).^m    
        for k=kr
            A[k,m] += C*(μ+k-1)
        end
    else
        C=2.^m.*factorial(μ-1)./factorial(λ-1)*(2./length(d)).^m        
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


rangespace{λ}(D::Integral{UltrasphericalSpace{λ}})=UltrasphericalSpace{λ-D.order}(domain(D))

function addentries!{λ}(D::Integral{UltrasphericalSpace{λ}},A::ShiftArray,kr::Range1)
    @assert λ==1
    @assert D.order==1   
    
    d=domain(D)
    
    for k=max(kr[1],2):kr[end]
        A[k,-1] += .5length(d)./(k-1)
    end
    
    A
end
 



## Conversion Operator




function Conversion{a,b}(A::UltrasphericalSpace{a},B::UltrasphericalSpace{b})
    @assert b > a

    if b==a+1
        Conversion{UltrasphericalSpace{a},UltrasphericalSpace{b},Float64}(A,B)
    else
        d=domain(A)
        Conversion(UltrasphericalSpace{b-1}(d),B)*Conversion(A,UltrasphericalSpace{b-1}(d))
    end
end   


function addentries!(M::Conversion{ChebyshevSpace,UltrasphericalSpace{1}},A::ShiftArray,kr::Range)
    for k=kr
        A[k,0] += (k == 1)? 1. : .5
        A[k,2] += -.5        
    end
    
    A    
end

function addentries!{m,λ}(M::Conversion{UltrasphericalSpace{m},UltrasphericalSpace{λ}},A::ShiftArray,kr::Range)
    @assert λ==m+1
    for k=kr
        A[k,0] += (λ-1.)./(k - 2. + λ)
        A[k,2] += -(λ-1.)./(k + λ)
    end
    
    A    
end

function multiplyentries!(M::Conversion{ChebyshevSpace,UltrasphericalSpace{1}},A::ShiftArray,kr::Range)
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

function multiplyentries!{m,λ}(M::Conversion{UltrasphericalSpace{m},UltrasphericalSpace{λ}},A::ShiftArray,kr::Range)
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

bandinds{m,λ}(C::Conversion{UltrasphericalSpace{m},UltrasphericalSpace{λ}})=0,2



## spaceconversion

# return the space that has banded Conversion to the other
function conversion_rule{aorder,border}(a::UltrasphericalSpace{aorder},b::UltrasphericalSpace{border})
    @assert domainscompatible(a,b)
    
    aorder < border?a:b
end


spaceconversion(g::Vector,::UltrasphericalSpace{1},::ChebyshevSpace)=ultraiconversion(g)
spaceconversion(g::Vector,::ChebyshevSpace,::UltrasphericalSpace{1})=ultraconversion(g)



## Dirichlet Conversion

addentries!(C::Conversion{ChebyshevDirichletSpace{1,0},ChebyshevSpace},A::ShiftArray,kr::Range1)=toeplitz_addentries!(ShiftVector([1.,1.],1),A,kr)
addentries!(C::Conversion{ChebyshevDirichletSpace{0,1},ChebyshevSpace},A::ShiftArray,kr::Range1)=toeplitz_addentries!(ShiftVector([1.,-1.],1),A,kr)
addentries!(C::Conversion{ChebyshevDirichletSpace{1,1},ChebyshevSpace},A::ShiftArray,kr::Range1)=toeplitz_addentries!(ShiftVector([1.,0.,-1.],1),A,kr)    
function addentries!(C::Conversion{ChebyshevDirichletSpace{2,2},ChebyshevSpace},A::ShiftArray,kr::Range1)
    for k=kr
        A[k,0]=1
        A[k,4]=2*(k+1)/k-1
        if k>= 3
            A[k,2]=-2*(k-1)/(k-2)
        end
    end
    
    A
end
bandinds(::Conversion{ChebyshevDirichletSpace{1,0},ChebyshevSpace})=0,1
bandinds(::Conversion{ChebyshevDirichletSpace{0,1},ChebyshevSpace})=0,1
bandinds(::Conversion{ChebyshevDirichletSpace{1,1},ChebyshevSpace})=0,2

conversion_rule(b::ChebyshevDirichletSpace,a::ChebyshevSpace)=b



## Evaluation Functional


function getindex(B::Evaluation{ChebyshevDirichletSpace{1,0},Bool},kr::Range)
    d = domain(B)
    
    if B.x == false && B.order == 0
        Float64[k==1?1.0:0.0 for k=kr]
    else
        getindex(Evaluation(d)*Conversion(domainspace(B)),kr)
    end
end

function getindex(B::Evaluation{ChebyshevDirichletSpace{0,1},Bool},kr::Range)
    d = domain(B)
    
    if B.x == true && B.order == 0
        Float64[k==1?1.0:0.0 for k=kr]
    else
        getindex(Evaluation(d)*Conversion(domainspace(B)),kr)
    end
end

function getindex(B::Evaluation{ChebyshevDirichletSpace{1,1},Bool},kr::Range)
   tol = 200.*eps()
    d = domain(B)
    
    if B.x == false && B.order == 0
        Float64[k==1?1.0:(k==2?-1.0:0.0) for k=kr]
    elseif B.x == true && B.order == 0
        Float64[k<=2?1.0:0.0 for k=kr]
    else
        getindex(Evaluation(d)*Conversion(domainspace(B)),kr)
    end
end
