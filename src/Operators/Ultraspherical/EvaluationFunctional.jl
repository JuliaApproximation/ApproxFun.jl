export dirichlet, neumann
export EvaluationFunctional, BasisFunctional

## EvaluationFunctional constructors

type EvaluationFunctional{D<:IntervalDomain,T<:Number} <: Functional{T}
    domain::D
    x::T
    order::Int
end

EvaluationFunctional(x::Number)=EvaluationFunctional(Interval(),x,0)
EvaluationFunctional(d::IntervalDomain,x::Number)=EvaluationFunctional(d,x,0)
EvaluationFunctional{T<:Number}(d::Vector{T},x::Number,o::Integer)=EvaluationFunctional(Interval(d),x,o)


domainspace(E::EvaluationFunctional)=ChebyshevSpace(E.domain)

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
function Base.getindex{D,T}(op::EvaluationFunctional{D,T},k::Range1)
   tol = 200.*eps()
    x = op.x
    d = op.domain

    ret = 0.
    
    if abs(x-d.a) < tol && op.order ==0
        ret=-(-1.).^k  ##TODO: speed up
    elseif  abs(x-d.b) < tol && op.order ==0
        ret=ones(size(k)[1])
    elseif  abs(x-d.a) < tol && op.order ==1
        ret=(k-1).*(k-1).*(-1.).^k*2/(d.b-d.a) 
    elseif  abs(x-d.b) < tol && op.order ==1
        ret=(k-1).*(k-1)*2/(d.b-d.a) 
    elseif  abs(x-d.a) < tol && op.order ==2
        ret=-(k.-1).^2.*((k.-1).^2.-1)/3.*(-1.).^k*(2/(d.b-d.a))^2
    elseif  abs(x-d.b) < tol && op.order ==2
        ret=(k.-1).^2.*((k.-1).^2.-1)/3*(2/(d.b-d.a))^2  
    elseif op.order == 0
        
        ret=evaluatechebyshev(k[end],tocanonical(d,x))[k]
    else
        error("Only zero–second order implemented")
    end
    
    ret
end


