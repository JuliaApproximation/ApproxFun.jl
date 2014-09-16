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

 




include("EvaluationFunctional.jl")
include("USConversionOperator.jl")
include("IntegrationOperator.jl")
include("DirichletConversionOperator.jl")


