export EvaluationOperator,DifferentialOperator,Operator
export differentialorder,conversionmatrix,derivativematrix,tomatrix
export toeplitzmatrix,hankelmatrix,multiplicationmatrix

abstract Operator


##Operator types

type EvaluationOperator <: Operator
    x
    domain::IntervalDomain
end

EvaluationOperator(x)=EvaluationOperator(x,Interval())
EvaluationOperator(x,d::Vector)=EvaluationOperator(x,Interval(d))

type DifferentialOperator <: Operator
    coefficients::Vector
    domain::IntervalDomain
end

#TODO ensure any funs match coefficients

DifferentialOperator(cfs)=DifferentialOperator(cfs,Interval())
DifferentialOperator(cfs,d::Vector)=DifferentialOperator(cfs,Interval(d))





##Information about operator


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
        p = zeros(k)
        p[1] = 1.
        p[2] = op.x
        
        for j=2:k-1
            p[j+1] = 2x*p[j] - p[j-1]
        end
        
        p[k]
  end
end

function Base.getindex(op::EvaluationOperator,k::Range1)
   tol = 200.*eps()

    if abs(op.x-op.domain.a) < tol
        -(-1.).^k
    elseif  abs(op.x-op.domain.b) < tol
        ones(size(k)[1])
    else
        #TODO k[end] == 1
        
        x = tocanonical(op.domain,op.x)
        
        p = zeros(k[end])
        p[1] = 1.
        p[2] = x
        
        for j=2:k[end]-1
            p[j+1] = 2x*p[j] - p[j-1]
        end
        
        p[k]
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
    m = max(k[end],j[end])  #TODO: smart indexing
    
    ret = spzeros(m,m)
    
    od = differentialorder(op)
    
    for l = 1:length(op.coefficients)
        ret += conversionmatrix(l-1:od,m)*multiplicationmatrix(l-1,op.coefficients[l],m)*derivativematrix(l-1,op.domain,m)
    end
    
    ret[k,j]
end



Base.getindex(op::DifferentialOperator,k::Integer,j::Integer)=op[k:k,j:j][1,1]

## Multiplication matrix

function toeplitzmatrix(v::Vector,n::Integer)
    ret = spzeros(n,n)

    for k=1:n
        ret[k,k] = 2v[1]
    end
    
    for j=1:length(v)-1
        for k=1:n-j
            ret[k,k + j] = v[j+1]
            ret[k + j,k] = v[j+1]            
        end
    end
  
    ret
end

function hankelmatrix(v::Vector,n::Integer)
    ret = spzeros(n,n)
    
    for j=1:length(v)
        for k=max(j-n+1,1):min(j,n)
            ret[k,j-k+1] = v[j]        
        end
    end
  
    ret
end


toeplitzmatrix(f::IFun,n::Integer)=toeplitzmatrix(f.coefficients,n)

multiplicationmatrix(f::IFun,n::Integer)=.5(toeplitzmatrix(f,n) + [spzeros(1,n),hankelmatrix(f.coefficients[2:end],n)[1:end-1,:]])

umultiplicationmatrix(f::IFun,n::Integer)=.5(toeplitzmatrix(f,n)  - hankelmatrix(f.coefficients[3:end],n))


multiplicationmatrix(a::Number,n::Integer)=a*speye(n)

function multiplicationmatrix(l::Integer,a,n)
    if typeof(a) <: Number
        return a*speye(n)
    elseif l==0
        multiplicationmatrix(a,n)
    elseif l==1
        umultiplicationmatrix(a,n)
    else
        error("higher order multiplication not yet implemented")
    end
end

## ultraconv

ultraconversion(v::Vector{Float64},μ::Integer)=conversionmatrix(0:μ,length(v))*v
ultraiconversion(v::Vector{Float64},μ::Integer)=conversionmatrix(0:μ,length(v))\v
coefficients(f::IFun,m::Integer)=ultraconversion(f.coefficients,m)

## Multiplication of operator * fun

#TODO: bet bottom length right
function *(A::DifferentialOperator,b::IFun)
    n=length(b)
    od=differentialorder(A)
    IFun(ultraiconversion(A[1:n,1:n]*b.coefficients,od),b.domain)
end


*(A::EvaluationOperator,b::IFun)=dot(A[1:length(b)],b.coefficients)
*(A::Array{Operator,1},b::IFun)=map(a->a*b,convert(Array{Any,1},A))

## Addidition of Differential operators

function +(A::DifferentialOperator,B::DifferentialOperator)
    @assert A.domain == B.domain
    
    n = max(length(A.coefficients),length(B.coefficients))
    
    DifferentialOperator(pad(A.coefficients,n) + pad(B.coefficients,n),A.domain)
end

-(A::DifferentialOperator)=DifferentialOperator(-A.coefficients,A.domain)

function -(A::DifferentialOperator,B::DifferentialOperator)
    @assert A.domain == B.domain
    
    n = max(length(A.coefficients),length(B.coefficients))
    
    DifferentialOperator(pad(A.coefficients,n) - pad(B.coefficients,n),A.domain)
end

*(a::Number,B::DifferentialOperator)=DifferentialOperator(a.*B.coefficients,B.domain)



## Linear Solve

function tomatrix(A::Array{Operator,1},n::Integer)  
  B = spzeros(n,n)
  
  nbc = mapreduce(a-> size(a)[1] == Inf ? 0 : size(a)[1],+,A)
  
  mapreduce(a-> size(a)[1] == Inf ? a[1:n-nbc,1:n] : a[1:size(a)[1],1:n],vcat,A)
end

\(A::Array{Operator,1},b::Array{Any,1})=\(A,b,eps())

function \(A::Array{Operator,1},b::Array{Any,1},tol::Float64)
    @assert length(A) == length(b)

    d=A[1].domain
    
    map(x -> (@assert d == x.domain),A)
    
    map(x -> typeof(x) <: IFun ? (@assert d == x.domain) : 0,b)  
    
    #min length
    m = mapreduce(x-> typeof(x) <: IFun ? length(x) : 1, +, b)
    
    #number of bc rows
    nbc=mapreduce(a-> size(a)[1] == Inf ? 0 : size(a)[1],+,A)    
    
    rhs=copy(b)
    
    ##TODO: relative vs absolute accuracy
#    nrm=mapreduce(norm,max,b)
    
    for logn = 3:20
        n = 2^logn + m
    
        M = tomatrix(A,n)
        

        #size of each part
        sz = map(a-> size(a)[1]== Inf ? n - nbc : size(a)[1],A)
        

        
        for j = 1:length(b)
            if typeof(b[j]) <: IFun
                rhs[j] = pad(coefficients(b[j],differentialorder(A)),sz[j]) 
            else
                rhs[j] = b[j].*ones(sz[j])
            end
        end

        
        cfs = M\vcat(rhs...)
        
        if (maximum(abs(cfs[end-8:end])) < tol)
            return IFun(chop(cfs,eps()),d)
        end
    end
end

\(A::Array{Operator,1},b::Union(Array{Float64,1},Array{Float64,2}))=A\convert(Array{Any,1},b)






