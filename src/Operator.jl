export EvaluationOperator,DifferentialOperator,Operator
export differentialorder,conversionmatrix,derivativematrix,tomatrix

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
Base.diff(d::IntervalDomain,n::Integer)=DifferentialOperator([zeros(n),1.],d)
Base.diff(d::IntervalDomain)=Base.diff(d,1)




function derivativematrix(μ::Integer,d::Domain,k::Range1,j::Range1)
    if(μ == 0)
        @assert k[1]==1 && j[1]==1
    
        speye(k[end],j[end])
    else
        ret = spzeros(size(k)[1],size(j)[1])
        
        st = max(k[1],j[1]-μ)
        en = min(k[end],j[end]-μ)
        
        for l=st:en
            ret[l - k[1] + 1, l - j[1] + 1 + μ] =  2.^(μ-1).*factorial(μ-1).*(2./(d.b-d.a)).^μ.*(μ+l-1)
        end
        
        ret
    end
end

function conversionmatrix(λ::Integer,k::Range1,j::Range1)
    ret = spzeros(size(k)[1],size(j)[1])

    if λ == 0
        for l=max(k[1],j[1]):min(k[end],j[end])
            ret[l - k[1] + 1, l - j[1] + 1] =  l==1? 1. : .5
        end
        
        for l=max(k[1],j[1]):min(k[end],j[end])-2
            ret[l - k[1] + 1, l - j[1] + 3] =  -.5
        end
        
        ret    
    else
        for l=max(k[1],j[1]):min(k[end],j[end])
            ret[l - k[1] + 1, l - j[1] + 1] =  λ./(l - k[1] + λ)
        end
        
        for l=max(k[1],j[1]):min(k[end],j[end])-2
            ret[l - k[1] + 1, l - j[1] + 3] =  -λ./(l - k[1] + λ + 2)
        end
        
        ret    
    end
end

function conversionmatrix(od::Range1,k::Range1,j::Range1)
    @assert k[1] == 1 && j[1] == 1
    
    m = max(k[end],j[end])
    
    ret = speye(m,m)
    
    for λ = od[1]:od[end]-1
       ret = conversionmatrix(λ,1:m,1:m)*ret 
    end
    
    ret[k,j]
end

function derivativematrix(μ::Integer,od::Range1,d::Domain,k::Range1,j::Range1)
    @assert od[1] == 0  #TODO: more general
    @assert k[1] == 1 && j[1] == 1
    
    conversionmatrix(μ:od[end],k,j)*derivativematrix(μ,d,k,j)
end


function Base.getindex(op::DifferentialOperator,k::Range1,j::Range1)
    @assert k[1] == 1 && j[1] == 1
    
    ret = spzeros(k[end],j[end])
    
    od = differentialorder(op)
    
    for l = 1:length(op.coefficients)
        ret += op.coefficients[l]*derivativematrix(l-1,0:od,op.domain,k,j)
    end
    
    ret
end



Base.getindex(op::DifferentialOperator,k::Integer,j::Integer)=op[k:k,j:j][1,1]


## ultraconv

ultraconv(v::Vector{Float64},μ::Integer)=conversionmatrix(0:μ,1:length(v),1:length(v))*v
coefficients(f::IFun,m::Integer)=ultraconv(f.coefficients,m)


## Linear Solve

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
