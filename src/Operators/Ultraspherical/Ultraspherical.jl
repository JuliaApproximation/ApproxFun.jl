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
function Base.getindex(op::EvaluationFunctional{Float64,Bool,ChebyshevSpace},k::Range)
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

function Base.getindex(op::EvaluationFunctional{Float64,Float64,ChebyshevSpace},k::Range)
    if op.order == 0    
        evaluatechebyshev(k[end],tocanonical(domain(op),op.x))[k]
    else
        error("Only zero–second order implemented")
    end
end



## Multiplication


function addentries!{T,D}(M::MultiplicationOperator{T,D,ChebyshevSpace},A::ShiftArray,kr::Range1)
    cfs=coefficients(M.f)
    toeplitz_addentries!(.5cfs,A,kr)
    hankel_addentries!(.5cfs,A,max(kr[1],2):kr[end])            
end

function addentries!{T,D}(M::MultiplicationOperator{T,D,UltrasphericalSpace{1}},A::ShiftArray,kr::Range1)
    cfs=coefficients(M.f)
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

function addentries!{T,D,λ}(M::MultiplicationOperator{T,D,UltrasphericalSpace{λ}},A::ShiftArray,kr::Range)
    a=coefficients(M.f,domainspace(M))
    for k=kr
        A[k,0]=a[1] 
    end

    if length(a) > 1
        jkr=max(1,kr[1]-length(a)+1):kr[end]+length(a)-1
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


#DerivativeOperator(k::Integer,d::IntervalDomain)=DerivativeOperator(k-1:k,d)
#DerivativeOperator(d::IntervalDomain)=DerivativeOperator(1,d)


rangespace{λ}(D::DerivativeOperator{Float64,UltrasphericalSpace{λ}})=UltrasphericalSpace{λ+D.order}(domain(D))

function addentries!{λ}(D::DerivativeOperator{Float64,UltrasphericalSpace{λ}},A::ShiftArray,kr::Range1)
    m=D.order
    d=domain(D)
    μ=λ+m

    @assert isa(d,Interval)

    if λ == 0
        C=2.^(m-1).*factorial(μ-1)*(2./(d.b-d.a)).^m    
        for k=kr
            A[k,m] += C*(μ+k-1)
        end
    else
        C=2.^m.*factorial(μ-1)./factorial(λ-1)*(2./(d.b-d.a)).^m        
        for k=kr        
            A[k,m] += C
        end
    end
    
    A
end
 
## TODO: reimplement

# function DerivativeOperator(order::Range1,d::IntervalDomain)
#     @assert order[1] == 0 && order[end] <= 2  ##TODO other orders
#     
#     Mp = Fun(x->tocanonicalD(d,x),d)
#     
#     if order[end] == 1
#         Mp*USDerivativeOperator(0:1)
#     elseif order[end] == 2
#         (Mp.^2)*USDerivativeOperator(0:2) + diff(Mp)*USDerivativeOperator(0:1)
#     end
# end

 




include("USConversionOperator.jl")
include("IntegrationOperator.jl")
include("DirichletConversionOperator.jl")


