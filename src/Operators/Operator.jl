export Operator,RowOperator,InfiniteOperator
export differentialorder,tomatrix
export bandrange



abstract Operator
abstract RowOperator <: Operator
abstract InfiniteOperator <: Operator

abstract ShiftOperator <: Operator #For biinfinite operators
abstract InfiniteShiftOperator <: ShiftOperator
abstract RowShiftOperator <: ShiftOperator


include("EvaluationOperator.jl")
include("DifferentialOperator.jl")


differentialorder(A::Array{Operator,1})=mapreduce(differentialorder,max,A)

Base.size(::InfiniteOperator)=[Inf,Inf]
Base.size(::RowOperator)=Any[1,Inf] #use Any vector so the 1 doesn't become a float
Base.size(op::Operator,k::Integer)=size(op)[k]


Base.getindex(op::InfiniteOperator,k::Integer,j::Integer)=op[k:k,j:j][1,1]
Base.getindex(op::InfiniteOperator,k::Integer,j::Range1)=op[k:k,j][1,:]
Base.getindex(op::InfiniteOperator,k::Range1,j::Integer)=op[k,j:j][:,1]

## Multiplication of operator * fun

#TODO: get bottom length right
function *(A::InfiniteOperator,b::IFun)
    n=length(b)
    od=differentialorder(A)
    IFun(ultraiconversion(A[1:n,1:n]*b.coefficients,od),b.domain)
end

*(A::RowOperator,b::IFun)=dot(A[1:length(b)],b.coefficients)
*{T<:Operator}(A::Vector{T},b::IFun)=map(a->a*b,convert(Array{Any,1},A))



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




include("specialfunctions.jl")

