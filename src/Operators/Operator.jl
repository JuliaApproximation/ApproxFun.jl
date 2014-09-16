export Operator,Functional,InfiniteOperator
export bandrange, linsolve, periodic
export ldirichlet,rdirichlet,lneumann,rneumann





abstract Operator{T} #T is the entry type, Float64 or Complex{Float64}
abstract Functional{T} <: Operator{T}
abstract InfiniteOperator{T} <: Operator{T}   #Infinite Operators have + range
abstract BandedBelowOperator{T} <: InfiniteOperator{T}
abstract BandedOperator{T} <: BandedBelowOperator{T}





##TODO: Why do we need BandedOperator to check for row>0?

## We assume operators are T->T
rangespace(A::Operator)=AnySpace()
domainspace(A::Operator)=AnySpace()
rangespace(A::Functional)=ScalarSpace()
domain(A::Operator)=domain(domainspace(A))




Base.size(::InfiniteOperator)=[Inf,Inf]
Base.size(::Functional)=Any[1,Inf] #use Any vector so the 1 doesn't become a float
Base.size(op::Operator,k::Integer)=size(op)[k]


## geteindex

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




index(b::BandedBelowOperator)=1-bandinds(b)[1]  # index is the equivalent of BandedArray.index



## Construct operators


ShiftArray{T<:Number}(B::Operator{T},k::Range1,j::Range1)=addentries!(B,sazeros(T,k,j),k)
ShiftArray(B::Operator,k::Range1)=ShiftArray(B,k,bandrange(B))
BandedArray(B::Operator,k::Range1)=BandedArray(B,k,(k[1]+bandinds(B)[1]):(k[end]+bandinds(B)[end]))
BandedArray(B::Operator,k::Range1,cs)=BandedArray(ShiftArray(B,k,bandrange(B)),cs)

include("ShiftOperator.jl")
include("linsolve.jl")

include("SpaceOperator.jl")

include("ToeplitzOperator.jl")

include("ConstantOperator.jl")
include("ConversionOperator.jl")
include("MultiplicationOperator.jl")
include("DerivativeOperator.jl")

include("Ultraspherical/Ultraspherical.jl")

include("SavedOperator.jl")
include("AlmostBandedOperator.jl")
include("adaptiveqr.jl")


include("OperatorAlgebra.jl")

include("specialfunctions.jl")

include("TransposeOperator.jl")
include("StrideOperator.jl")
include("SliceOperator.jl")

include("CompactOperator.jl")

include("Fourier/FourierDerivativeOperator.jl")

include("null.jl")
include("systems.jl")



## Convenience routines

Base.diff(d::DomainSpace,μ::Integer)=DerivativeOperator(d,μ)
Base.diff(d::IntervalDomain,μ::Integer)=diff(ChebyshevSpace(d),μ)
Base.diff(d::PeriodicDomain,μ::Integer)=diff(LaurentSpace(d),μ)
Base.diff(d::Domain)=Base.diff(d,1)


Base.zero{T<:Number}(::Type{Operator{T}})=ConstantOperator(zero(T))
Base.zero{O<:Operator}(::Type{O})=ConstantOperator(0.0)

integrate(d::IntervalDomain)=IntegrationOperator(1,d)

evaluate(d::IntervalDomain,x)=EvaluationFunctional(d,x)
ldirichlet(d::IntervalDomain)=evaluate(d,d.a)
rdirichlet(d::IntervalDomain)=evaluate(d,d.b)
lneumann(d::IntervalDomain)=EvaluationFunctional(d,d.a,1)
rneumann(d::IntervalDomain)=EvaluationFunctional(d,d.b,1)


for func in (:ldirichlet,:rdirichlet,:lneumann,:rneumann)
    @eval ($func)(d::IntervalDomainSpace)=$func(domain(d))*ConversionOperator(d)
end


dirichlet(d::IntervalDomain)=[ldirichlet(d),rdirichlet(d)]
dirichlet(d::IntervalDomainSpace)=[ldirichlet(d),rdirichlet(d)]
neumann(d::IntervalDomain)=[lneumann(d),rneumann(d)]
neumann(d::IntervalDomainSpace)=[lneumann(d),rneumann(d)]


function dirichlet{T<:IntervalDomain}(d::Vector{T})
    m=length(d)
    B=zeros(Operator,2m,m)
    B[1,1]=dirichlet(d[1])[1]
    B[2,end]=dirichlet(d[end])[end]
    for k=1:m-1
        B[k+2,k]=dirichlet(d[k])[2]
        B[k+2,k+1]=-dirichlet(d[k+1])[1]    
        B[k+m+1,k]=neumann(d[k])[2]
        B[k+m+1,k+1]=-neumann(d[k+1])[1]        
    end
    B
end

function periodic{T<:IntervalDomain}(d::Vector{T})
    m=length(d)
    B=zeros(Operator,2m,m)
    B[1,1]=dirichlet(d[1])[1]
    B[1,end]=-dirichlet(d[end])[end]
    B[2,1]=neumann(d[1])[1]
    B[2,end]=-neumann(d[end])[end]

    for k=1:m-1
        B[k+2,k]=dirichlet(d[k])[2]
        B[k+2,k+1]=-dirichlet(d[k+1])[1]    
        B[k+m+1,k]=neumann(d[k])[2]
        B[k+m+1,k+1]=-neumann(d[k+1])[1]        
    end
    B
end

## Conversion

# TODO: can convert return different type?
Base.convert{T<:Operator}(A::Type{T},n::Number)=ConstantOperator(1.0*n)


## Promotion

for T in (:Float64,:Int64,:(Complex{Float64}))
    @eval Base.promote_rule{N<:Number,O<:Operator{$T}}(::Type{N},::Type{O})=Operator{promote_type(N,$T)}
end

Base.promote_rule{N<:Number,O<:Operator}(::Type{N},::Type{O})=Operator
