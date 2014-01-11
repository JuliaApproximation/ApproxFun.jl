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


ultraiconversion(g::Vector,m::Integer)=backsubstitution!(MutableAlmostBandedOperator(Operator[ConversionOperator(0:m)]),copy(g))
ultraconversion(g::Vector,m::Integer)=ConversionOperator(0:m)*g

*(A::BandedBelowOperator,b::Vector)= A[1:length(b)-bandrange(A)[1],1:length(b)]*b
*(A::InfiniteOperator,b::IFun)=IFun(ultraiconversion(A*b.coefficients,rangespace(A)),b.domain)

*(A::RowOperator,b::Vector)=dot(A[1:length(b)],b)
*(A::RowOperator,b::IFun)=A*b.coefficients
*{T<:Operator}(A::Vector{T},b::IFun)=map(a->a*b,convert(Array{Any,1},A))



## Linear Solve


##TODO: reimplement tolerance
\(A::Vector{Operator},b::Vector)=IFun(adaptiveqr(A,b),A[1].domain)



include("ShiftArray.jl")


include("ToeplitzOperator.jl")
include("EvaluationOperator.jl")
include("ConversionOperator.jl")
include("DerivativeOperator.jl")

include("AlmostBandedOperator.jl")
include("OperatorAlgebra.jl")

include("specialfunctions.jl")

