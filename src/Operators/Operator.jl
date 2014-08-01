export Operator,Functional,InfiniteOperator
export bandrange, linsolve





abstract Operator{T} #T is the entry type, Flaot64 or Complex{Float64}
abstract Functional{T} <: Operator{T}
abstract InfiniteOperator{T} <: Operator{T}   #Infinite Operators have + range
abstract BandedBelowOperator{T} <: InfiniteOperator{T}
abstract BandedOperator{T} <: BandedBelowOperator{T}

abstract ShiftOperator{T} <: Operator{T} #For biinfinite operators
abstract InfiniteShiftOperator{T} <: ShiftOperator{T}
abstract BandedShiftOperator{T} <: InfiniteShiftOperator{T}
abstract ShiftFunctional{T} <: Functional{T}


##TODO: Change BandedOperator -> BandedInfiniteOperator
##TODO: Remove InfiniteShiftOperator
##TODO: Add BandedOperator = Union(BandedInfiniteOperator,BandedShiftOperator)
##TODO: Why do we need BandedOperator to check for row>0?

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

commondomain{T<:Number}(P::Vector,g::Array{T})=commondomain(P)
commondomain(P::Vector,g)=commondomain([P,g])


Base.size(::InfiniteOperator)=[Inf,Inf]
Base.size(::Functional)=Any[1,Inf] #use Any vector so the 1 doesn't become a float
Base.size(op::Operator,k::Integer)=size(op)[k]


Base.getindex(op::Operator,k::Integer,j::Integer)=op[k:k,j:j][1,1]
Base.getindex(op::Operator,k::Integer,j::Range1)=op[k:k,j][1,:]
Base.getindex(op::Operator,k::Range1,j::Integer)=op[k,j:j][:,1]


Base.getindex(op::Functional,k::Integer)=op[k:k][1]

function Base.getindex(op::Functional,j::Range1,k::Range1)
  @assert j[1]==1 && j[end]==1
  op[k]' #TODO conjugate transpose?
end
function Base.getindex(op::Functional,j::Integer,k::Range1)
  @assert j==1
  op[k]' #TODO conjugate transpose?
end



function Base.getindex(B::Operator,k::Range1,j::Range1)
    BandedArray(B,k,j)[k,j]
end


## bandrange and indexrange

bandrange(b::BandedShiftOperator)=Range1(bandinds(b)...)
bandrange(b::BandedBelowOperator)=Range1(bandinds(b)...)
function bandrangelength(B::BandedBelowOperator)
    bndinds=bandinds(B)
    bndinds[end]-bndinds[1]+1
end


function columninds(b::BandedBelowOperator,k::Integer)
    ret = bandinds(b)
  
    (ret[1]  + k < 1) ? (1,(ret[end] + k)) : (ret[1]+k,ret[2]+k)
end

##TODO: Change to columnindexrange to match BandedOperator
indexrange(b::BandedBelowOperator,k::Integer)=Range1(columninds(b,k)...)




index(b::BandedBelowOperator)=1-bandinds(b)[1]



## Construct operators


ShiftArray{T<:Number}(B::Operator{T},k::Range1,j::Range1)=addentries!(B,sazeros(T,k,j),k)
ShiftArray(B::Operator,k::Range1)=ShiftArray(B,k,bandrange(B))
BandedArray(B::Operator,k::Range1)=BandedArray(B,k,(k[1]+bandinds(B)[1]):(k[end]+bandinds(B)[end]))
BandedArray(B::Operator,k::Range1,cs)=BandedArray(ShiftArray(B,k,bandrange(B)),cs)


include("linsolve.jl")

include("OperatorSpace.jl")

include("ToeplitzOperator.jl")

include("ConstantOperator.jl")

include("Ultraspherical/Ultraspherical.jl")

include("SavedOperator.jl")
include("AlmostBandedOperator.jl")
include("adaptiveqr.jl")


include("OperatorAlgebra.jl")
include("FunctionalAlgebra.jl")

include("specialfunctions.jl")

include("TransposeOperator.jl")
include("StrideOperator.jl")
include("SliceOperator.jl")

include("CompactOperator.jl")

include("Fourier/FourierDerivativeOperator.jl")
include("Fourier/FourierSpace.jl")

include("null.jl")
include("systems.jl")



## Convenience routines

Base.diff(d::IntervalDomain,μ::Integer)=DerivativeOperator(0:μ,d)
Base.diff(d::PeriodicDomain,μ::Integer)=FourierDerivativeOperator(μ,d)
Base.diff(d::Domain)=Base.diff(d,1)

Base.eye(d::IntervalDomain)=MultiplicationOperator(IFun([1.],d))
Base.eye(d::PeriodicDomain)=MultiplicationOperator(FFun(ShiftVector([1.],1),d))
Base.zero{T<:Number}(::Type{Operator{T}})=ConstantOperator(zero(T))
Base.zero{O<:Operator}(::Type{O})=ConstantOperator(0.0)

integrate(d::IntervalDomain)=IntegrationOperator(1,d)

evaluate(d::IntervalDomain,x)=EvaluationFunctional(d,x)
dirichlet(d::IntervalDomain)=[evaluate(d,d.a),evaluate(d,d.b)]
neumann(d::IntervalDomain)=[EvaluationFunctional(d,d.a,1),EvaluationFunctional(d,d.b,1)]


## Conversion

Base.convert{T<:Number}(A::Type{Operator{T}},n::Number)=ConstantOperator(one(T)*n)

## Promotion

for T in (:Float64,:Int64,:(Complex{Float64}))
    @eval Base.promote_rule{N<:Number,O<:Operator{$T}}(::Type{N},::Type{O})=Operator{promote_type(N,$T)}
end

Base.promote_rule{N<:Number,O<:Operator}(::Type{N},::Type{O})=Operator
