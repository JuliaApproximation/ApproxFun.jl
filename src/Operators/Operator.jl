export Operator,RowOperator,InfiniteOperator
export bandrange



abstract Operator
abstract RowOperator <: Operator
abstract InfiniteOperator <: Operator
abstract BandedBelowOperator <: InfiniteOperator
abstract BandedOperator <: BandedBelowOperator

abstract ShiftOperator <: Operator #For biinfinite operators
abstract InfiniteShiftOperator <: ShiftOperator
abstract RowShiftOperator <: ShiftOperator


## We assume operators are T->T
rangespace(A::InfiniteOperator)=0
domainspace(A::InfiniteOperator)=0
domain(A::Operator)=Any

domain(f::IFun)=f.domain
domain(::Number)=Any

function domain(P::Vector)
    ret = Any
    
    for op in P
        d = domain(op)
        @assert ret == Any || d == Any || ret == d
        
        if d != Any
            ret = d
        end
    end
    
    ret
end



Base.size(::InfiniteOperator)=[Inf,Inf]
Base.size(::RowOperator)=Any[1,Inf] #use Any vector so the 1 doesn't become a float
Base.size(op::Operator,k::Integer)=size(op)[k]


Base.getindex(op::InfiniteOperator,k::Integer,j::Integer)=op[k:k,j:j][1,1]
Base.getindex(op::InfiniteOperator,k::Integer,j::Range1)=op[k:k,j][1,:]
Base.getindex(op::InfiniteOperator,k::Range1,j::Integer)=op[k,j:j][:,1]


function Base.getindex(B::BandedOperator,k::Range1,j::Range1)
    BandedArray(B,k)[k,j]
end



## Multiplication of operator * fun


ultraiconversion(g::Vector,m::Integer)=(m==0)? g : backsubstitution!(MutableAlmostBandedOperator(Operator[ConversionOperator(0:m)]),copy(g))
ultraconversion(g::Vector,m::Integer)=(m==0)? g : ConversionOperator(0:m)*g

*(A::BandedBelowOperator,b::Vector)= A[1:length(b)-bandrange(A)[1],1:length(b)]*b
*(A::InfiniteOperator,b::IFun)=IFun(ultraiconversion(A*ultraconversion(b.coefficients,domainspace(A)),rangespace(A)),b.domain)

*(A::RowOperator,b::Vector)=dot(A[1:length(b)],b)
*(A::RowOperator,b::IFun)=A*b.coefficients
*{T<:Operator}(A::Vector{T},b::IFun)=map(a->a*b,convert(Array{Any,1},A))



## Linear Solve



\{T<:Operator}(A::Vector{T},b::Vector,tol::Float64)=IFun(adaptiveqr(A,b,tol),domain([A,b]))
\{T<:Operator}(A::Vector{T},b::Vector)=IFun(adaptiveqr(A,b),domain([A,b]))
\(A::Operator,b::Vector)=[A]\b
\(A::Operator,b::IFun)=[A]\[b]



include("ShiftArray.jl")


include("ToeplitzOperator.jl")
include("EvaluationOperator.jl")
include("ConversionOperator.jl")
include("DerivativeOperator.jl")
include("IntegrationOperator.jl")

include("AlmostBandedOperator.jl")
include("OperatorAlgebra.jl")

include("specialfunctions.jl")

