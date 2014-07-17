export Operator,RowOperator,InfiniteOperator
export bandrange, linsolve



abstract Operator{T} #T is the entry type, Flaot64 or Complex{Float64}
abstract RowOperator{T} <: Operator{T}
abstract InfiniteOperator{T} <: Operator{T}
abstract BandedBelowOperator{T} <: InfiniteOperator{T}
abstract BandedOperator{T} <: BandedBelowOperator{T}

abstract ShiftOperator{T} <: Operator{T} #For biinfinite operators
abstract InfiniteShiftOperator{T} <: ShiftOperator{T}
abstract BandedShiftOperator{T} <: InfiniteShiftOperator{T}
abstract RowShiftOperator{T} <: ShiftOperator{T}


## We assume operators are T->T
domain(A::Operator)=Any


function commondomain(P::Vector)
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

commondomain{T<:Number}(P::Vector,g::Vector{T})=commondomain(P)
commondomain(P::Vector,g)=commondomain([P,g])


Base.size(::InfiniteOperator)=[Inf,Inf]
Base.size(::RowOperator)=Any[1,Inf] #use Any vector so the 1 doesn't become a float
Base.size(op::Operator,k::Integer)=size(op)[k]


Base.getindex(op::Operator,k::Integer,j::Integer)=op[k:k,j:j][1,1]
Base.getindex(op::Operator,k::Integer,j::Range1)=op[k:k,j][1,:]
Base.getindex(op::Operator,k::Range1,j::Integer)=op[k,j:j][:,1]


Base.getindex(op::RowOperator,k::Integer)=op[k:k][1]

function Base.getindex(op::RowOperator,j::Range1,k::Range1)
  @assert j[1]==1 && j[end]==1
  op[k]' #TODO conjugate transpose?
end
function Base.getindex(op::RowOperator,j::Integer,k::Range1)
  @assert j==1
  op[k]' #TODO conjugate transpose?
end



function Base.getindex(B::Operator,k::Range1,j::Range1)
    BandedArray(B,k,j)[k,j]
end


## indexrange

function indexrange(b::BandedBelowOperator,k::Integer)
    ret = bandrange(b) + k
  
    (ret[1] < 1) ? (1:ret[end]) : ret
end

index(b::BandedBelowOperator)=1-bandrange(b)[1]



## Construct operators


ShiftArray{T<:Number}(B::Operator{T},k::Range1,j::Range1)=addentries!(B,sazeros(T,k,j),k)
ShiftArray(B::Operator,k::Range1)=ShiftArray(B,k,bandrange(B))
BandedArray(B::Operator,k::Range1)=BandedArray(B,k,(k[1]+bandrange(B)[1]):(k[end]+bandrange(B)[end]))
BandedArray(B::Operator,k::Range1,cs)=BandedArray(ShiftArray(B,k,bandrange(B)),cs)


include("linsolve.jl")

include("OperatorSpace.jl")

include("ToeplitzOperator.jl")

include("ConstantOperator.jl")

include("Ultraspherical/Ultraspherical.jl")

include("SavedBandedOperator.jl")
include("AlmostBandedOperator.jl")
include("adaptiveqr.jl")


include("OperatorAlgebra.jl")
include("RowOperatorAlgebra.jl")

include("specialfunctions.jl")

include("TransposeOperator.jl")
include("StrideOperator.jl")
include("SliceOperator.jl")

include("CompactOperator.jl")

include("Fourier/FourierDerivativeOperator.jl")
include("Fourier/FourierSpace.jl")

include("null.jl")



## Convenience routines

Base.diff(d::IntervalDomain,μ::Integer)=DerivativeOperator(0:μ,d)
Base.diff(d::PeriodicDomain,μ::Integer)=FourierDerivativeOperator(μ,d)
Base.diff(d::Domain)=Base.diff(d,1)

Base.eye(d::IntervalDomain)=MultiplicationOperator(IFun([1.],d))
Base.eye(d::PeriodicDomain)=MultiplicationOperator(FFun(ShiftVector([1.],1),d))

integrate(d::IntervalDomain)=IntegrationOperator(1,d)

evaluate(d::IntervalDomain,x)=EvaluationOperator(d,x)
dirichlet(d::IntervalDomain)=[evaluate(d,d.a),evaluate(d,d.b)]
neumann(d::IntervalDomain)=[EvaluationOperator(d,d.a,1),EvaluationOperator(d,d.b,1)]

Base.start(d::IntervalDomain)=evaluate(d,d.a)
Base.endof(d::IntervalDomain)=evaluate(d,d.b)

## Conversion

Base.convert{N<:Number}(A::Type{Operator},n::N)=ConstantOperator(1.n)
