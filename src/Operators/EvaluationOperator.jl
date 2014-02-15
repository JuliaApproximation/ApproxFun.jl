export dirichlet, neumann
export EvaluationOperator, BasisOperator

## EvaluationOperator constructors

type EvaluationOperator{D<:IntervalDomain,T<:Number} <: RowOperator{T}
    domain::D
    x::T
    order::Integer
end

EvaluationOperator(x)=EvaluationOperator(Interval(),x,0)
EvaluationOperator(d,x)=EvaluationOperator(d,x,0)
EvaluationOperator(d::Vector,x,o)=EvaluationOperator(Interval(d),x,o)




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
function Base.getindex{D,T}(op::EvaluationOperator{D,T},k::Range1)
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
    elseif op.order == 0
        
        ret=evaluatechebyshev(k[end],tocanonical(d,x))[k]
    else
        error("Only first and zero order implemented")
    end
    
    ret
end


function Base.getindex(op::EvaluationOperator,j::Range1,k::Range1)
  @assert j[1]==1 && j[end]==1
  op[k]' #TODO conjugate transpose?
end
function Base.getindex(op::EvaluationOperator,j::Integer,k::Range1)
  @assert j==1
  op[k]' #TODO conjugate transpose?
end




type BasisOperator <: RowOperator{Float64}
    k::Integer
end


##TODO: the overloading as both vector and row vector may be confusing
Base.getindex(op::BasisOperator,k::Integer)=(k==op.k)?1.:0.

function Base.getindex(op::BasisOperator,k::Range1)
    convert(Vector{Float64},k.==op.k)
end


function Base.getindex(op::BasisOperator,j::Range1,k::Range1)
  @assert j[1]==1 && j[end]==1
  op[k]' #TODO conjugate transpose?
end
function Base.getindex(op::BasisOperator,j::Integer,k::Range1)
  @assert j==1
  op[k]' #TODO conjugate transpose?
end