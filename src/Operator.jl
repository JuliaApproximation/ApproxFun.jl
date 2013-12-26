export EvaluationOperator,DifferentialOperator,Operator

abstract Operator


type EvaluationOperator <: Operator
    x
    domain::IntervalDomain
end

EvaluationOperator(x)=EvaluationOperator(x,Interval())


type DifferentialOperator <: Operator
    coefficients::Vector
    domain::IntervalDomain
end


Base.size(::EvaluationOperator)=Any[1,Inf]
Base.size(::DifferentialOperator)=[Inf,Inf]


differentialorder(::EvaluationOperator)=0
differentialorder(d::DifferentialOperator)=length(d.coefficients)-1
differentialorder(A::Array{Operator,1})=mapreduce(differentialorder,max,A)

function Base.getindex(op::EvaluationOperator,k::Integer)
    tol = 200.*eps()

    if abs(op.x-op.domain.a) < tol
        (-1.)^(k-1)
    elseif  abs(op.x-op.domain.b) < tol
        1.
    else
        error("getindex not implemented for general")
  end
end

Base.getindex(op::EvaluationOperator,k::Range1)=[op[l] for l=k]
function Base.getindex(op::EvaluationOperator,j::Range1,k::Range1)
  @assert j[1]==1 && j[end]==1
  op[k]'
end


MultiplicationOperator(a::IFun)=DifferentialOperator([a])
Base.diff(d::IntervalDomain)=DifferentialOperator([1.,0],d)


function Base.getindex(op::DifferentialOperator,k::Integer,j::Integer)
  if (k+1==j)
    k
  else
    0
  end
end

function Base.getindex(op::DifferentialOperator,k::Range1,j::Range1)
  ret = spzeros(size(k)[1],size(j)[1])
  
  st = max(k[1],j[1]-1)
  en = min(k[end],j[end]-1)
  
  for l=st:en
    ret[l - k[1] + 1, l - j[1] + 2] = l
  end
  
  ret
end

function tomatrix(A::Array{Operator,1},n::Integer)  
  B = spzeros(n,n)
  
  nbc = mapreduce(a-> size(a)[1] == Inf ? 0 : size(a)[1],+,A)
  
  mapreduce(a-> size(a)[1] == Inf ? a[1:n-nbc,1:n] : a[1:size(a)[1],1:n],vcat,A)
end

function \(A::Array{Operator,1},b::Array{Any,1})
    d=A[1].domain
    
    map(x -> (@assert d == x.domain),A)
    
    map(x -> typeof(x) <: IFun ? (@assert d == x.domain) : 0,b)  
    
    
    m = mapreduce(x-> typeof(x) <: IFun ? length(x) : 1, +, b)
    
    for logn = 1:20
        n = 2^logn + m
    
        M = tomatrix(A,n)
        nbc=mapreduce(a-> size(a)[1] == Inf ? 0 : size(a)[1],+,A)
        
        
        
        
        rhs=mapreduce(x->typeof(x) <: IFun ? pad(coefficients(x,differentialorder(A)),n-nbc) : [x], vcat, b)
        
        cfs = M\rhs
        
        if (maximum(abs(cfs[end-8:end])) < eps())
            return IFun(chop(cfs,eps()),d)
        end
    end
end
