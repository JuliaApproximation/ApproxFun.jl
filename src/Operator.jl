export EvaluationOperator,DifferentialOperator,Operator
export differentialorder,conversionmatrix,derivativematrix,tomatrix
export toeplitzmatrix,hankelmatrix,multiplicationmatrix
export dirichlet, neumann

abstract Operator


## EvaluationOperator constructors

type EvaluationOperator <: Operator
    x
    domain::IntervalDomain
    order::Integer
end

EvaluationOperator(x)=EvaluationOperator(x,Interval(),0)
EvaluationOperator(x,d)=EvaluationOperator(x,d,0)
EvaluationOperator(x,d::Vector,o)=EvaluationOperator(x,Interval(d),o)


evaluate(d::IntervalDomain,x)=EvaluationOperator(x,d)
dirichlet(d::IntervalDomain)=[evaluate(d,d.a),evaluate(d,d.b)]
neumann(d::IntervalDomain)=[EvaluationOperator(d.a,d,1),EvaluationOperator(d.b,d,1)]



## DifferentialOperator constructors

type DifferentialOperator{T<:IFun} <: Operator
    coefficients::Vector{T}
    domain::IntervalDomain
end

##TODO: ensure any funs match coefficients

DifferentialOperator(cfs::Vector)=DifferentialOperator(cfs,Interval())
DifferentialOperator(cfs,d::Vector)=DifferentialOperator(cfs,Interval(d))
DifferentialOperator{T<:Real}(cfs::Vector{T},d::IntervalDomain)=DifferentialOperator(IFun([1.],d)*cfs,d)
DifferentialOperator(cfs::Vector{Any},d::IntervalDomain)=DifferentialOperator(IFun([1.],d)*cfs,d)


MultiplicationOperator(a::IFun)=DifferentialOperator([a])
Base.diff(d::IntervalDomain,n::Integer)=DifferentialOperator([zeros(n),1.],d)
Base.diff(d::IntervalDomain)=Base.diff(d,1)



##Information about operator


Base.size(::EvaluationOperator)=Any[1,Inf]
Base.size(::DifferentialOperator)=[Inf,Inf]


differentialorder(d::EvaluationOperator)=0 ##TODO: decide whether 0 is the correct
differentialorder(d::DifferentialOperator)=length(d.coefficients)-1
differentialorder(A::Array{Operator,1})=mapreduce(differentialorder,max,A)


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

multiplicationmatrix(l::Integer,a::Number,n)=a*speye(n)

function multiplicationmatrix(l::Integer,a::IFun,n)
    if length(a) == 1
        multiplicationmatrix(a.coefficients[1],n)
    elseif l==0
        multiplicationmatrix(a,n)
    elseif l==1
        umultiplicationmatrix(a,n)
    elseif 
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
*{T<:Operator}(A::Vector{T},b::IFun)=map(a->a*b,convert(Array{Any,1},A))

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

function +(a::IFun, B::DifferentialOperator)
    @assert a.domain == B.domain

    cfs = deepcopy(B.coefficients)
    cfs[1] = cfs[1] + a
    DifferentialOperator(cfs,B.domain)
end

+(B::DifferentialOperator,a::IFun)=a+B
-(B::DifferentialOperator,a::IFun)=-a + B


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






