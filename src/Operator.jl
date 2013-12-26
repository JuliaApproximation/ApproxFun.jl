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




function derivativematrix(μ::Integer,d::Domain,m::Integer)
    if(μ == 0)    
        speye(m,m)
    else
        ret = spzeros(m,m)

        for l=1:m-μ
            ret[l, l + μ] =  2.^(μ-1).*factorial(μ-1).*(2./(d.b-d.a)).^μ.*(μ+l-1)
        end
        
        ret
    end
end

function conversionmatrix(λ::Integer,m::Integer)
    ret = spzeros(m,m)

    if λ == 0
        for l=1:m
            ret[l, l] =  l==1? 1. : .5
        end
        
        for l=1:m-2
            ret[l, l + 2] =  -.5
        end
        
        ret    
    else
        for l=1:m
            ret[l, l] =  λ./(l - 1 + λ)
        end
        
        for l=1:m-2
            ret[l, l + 2] =  -λ./(l + λ + 1)
        end
        
        ret    
    end
end

function conversionmatrix(od::Range1,m::Integer)        
    ret = speye(m,m)
    
    for λ = od[1]:od[end]-1
       ret = conversionmatrix(λ,m)*ret 
    end
    
    ret
end

function derivativematrix(μ::Integer,od::Range1,d::Domain,m::Integer)
    @assert od[1] == 0  #TODO: more general
    
    conversionmatrix(μ:od[end],m)*derivativematrix(μ,d,m)
end


function Base.getindex(op::DifferentialOperator,k::Range1,j::Range1)
    m = max(k[end],j[end])
    
    ret = spzeros(m,m)
    
    od = differentialorder(op)
    
    for l = 1:length(op.coefficients)
        ret += op.coefficients[l]*derivativematrix(l-1,0:od,op.domain,m)
    end
    
    ret[k,j]
end



Base.getindex(op::DifferentialOperator,k::Integer,j::Integer)=op[k:k,j:j][1,1]


## ultraconv

ultraconv(v::Vector{Float64},μ::Integer)=conversionmatrix(0:μ,length(v))*v
coefficients(f::IFun,m::Integer)=ultraconv(f.coefficients,m)


## Linear Solve

function tomatrix(A::Array{Operator,1},n::Integer)  
  B = spzeros(n,n)
  
  nbc = mapreduce(a-> size(a)[1] == Inf ? 0 : size(a)[1],+,A)
  
  mapreduce(a-> size(a)[1] == Inf ? a[1:n-nbc,1:n] : a[1:size(a)[1],1:n],vcat,A)
end

function \(A::Array{Operator,1},b::Array{Any,1})
    @assert length(A) == length(b)

    d=A[1].domain
    
    map(x -> (@assert d == x.domain),A)
    
    map(x -> typeof(x) <: IFun ? (@assert d == x.domain) : 0,b)  
    
    
    m = mapreduce(x-> typeof(x) <: IFun ? length(x) : 1, +, b)
    
    rhs=copy(b)
    
    for logn = 3:20
        n = 2^logn + m
    
        M = tomatrix(A,n)
        
        #number of bc rows
        nbc=mapreduce(a-> size(a)[1] == Inf ? 0 : size(a)[1],+,A)
        #size of each part
        sz = map(a-> size(a)[1]== Inf ? n - nbc : size(a)[1],A)
        

        
        for j = 1:length(b)
            if typeof(b[j]) == IFun
                rhs[j] = pad(coefficients(x,differentialorder(A)),n-nbc) 
            else
                rhs[j] = b[j].*ones(sz[j])
            end
        end

        
        cfs = M\vcat(rhs...)
        
        if (maximum(abs(cfs[end-8:end])) < eps())
            return IFun(chop(cfs,eps()),d)
        end
    end
end
