export dirichlet, neumann
export EvaluationOperator

## EvaluationOperator constructors

type EvaluationOperator{D<:IntervalDomain,T<:Number} <: RowOperator
    domain::D
    x::T
    order::Integer
end

EvaluationOperator(x)=EvaluationOperator(Interval(),x,0)
EvaluationOperator(d,x)=EvaluationOperator(d,x,0)
EvaluationOperator(d::Vector,x,o)=EvaluationOperator(Interval(d),x,o)


evaluate(d::IntervalDomain,x)=EvaluationOperator(d,x)
dirichlet(d::IntervalDomain)=[evaluate(d,d.a),evaluate(d,d.b)]
neumann(d::IntervalDomain)=[EvaluationOperator(d,d.a,1),EvaluationOperator(d,d.b,1)]


differentialorder(d::EvaluationOperator)=0 ##TODO: decide whether 0 is the correct



function evaluatechebyshev(n::Integer,x)
    if n == 1
        [1.]
    else
        p = zeros(n)
        p[1] = 1.
        p[2] = x
        
        for j=2:n-1
            p[j+1] = 2x*p[j] - p[j-1]
        end
        
        p
    end
end


Base.getindex(op::EvaluationOperator,k::Integer)=op[k:k][1]


##TODO: the overloading as both vector and row vector may be confusing
function Base.getindex(op::EvaluationOperator,k::Range1)
   tol = 200.*eps()
    x = op.x
    d = op.domain

    if abs(x-d.a) < tol && op.order ==0
        -(-1.).^k
    elseif  abs(x-d.b) < tol && op.order ==0
        ones(size(k)[1])
    elseif  abs(x-d.a) < tol && op.order ==1
        (k-1).*(k-1).*(-1.).^k*2/(d.b-d.a) 
    elseif  abs(x-d.b) < tol && op.order ==1
        (k-1).*(k-1)*2/(d.b-d.a) 
    elseif op.order == 0
        
        evaluatechebyshev(k[end],tocanonical(d,x))[k]
    else
        error("Only first and zero order implemented")
    end
end


function Base.getindex(op::EvaluationOperator,j::Range1,k::Range1)
  @assert j[1]==1 && j[end]==1
  op[k]' #TODO conjugate transpose?
end
function Base.getindex(op::EvaluationOperator,j::Integer,k::Range1)
  @assert j==1
  op[k]' #TODO conjugate transpose?
end


