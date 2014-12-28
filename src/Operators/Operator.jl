export Operator,Functional,InfiniteOperator
export bandrange, linsolve, periodic
export dirichlet, neumann
export ldirichlet,rdirichlet,lneumann,rneumann
export ldiffbc,rdiffbc,diffbcs




abstract Operator{T} #T is the entry type, Float64 or Complex{Float64}
abstract Functional{T} <: Operator{T}
abstract InfiniteOperator{T} <: Operator{T}   #Infinite Operators have + range
abstract BandedBelowOperator{T} <: InfiniteOperator{T}
abstract BandedOperator{T} <: BandedBelowOperator{T}

Base.eltype{T}(::Operator{T})=T




## We assume operators are T->T
rangespace(A::Operator)=AnySpace()
domainspace(A::Operator)=AnySpace()
rangespace(A::Functional)=ConstantSpace()
domain(A::Operator)=domain(domainspace(A))




Base.size(::InfiniteOperator)=[Inf,Inf]
Base.size(::Functional)=Any[1,Inf] #use Any vector so the 1 doesn't become a float
Base.size(op::Operator,k::Integer)=size(op)[k]







## bandrange and indexrange


bandinds(A,k::Integer)=bandinds(A)[k]
bandrange(b::BandedBelowOperator)=Range1(bandinds(b)...)
function bandrangelength(B::BandedBelowOperator)
    bndinds=bandinds(B)
    bndinds[end]-bndinds[1]+1
end


function columninds(b::BandedBelowOperator,k::Integer)
    ret = bandinds(b)
  
    (ret[1]  + k < 1) ? (1,(ret[end] + k)) : (ret[1]+k,ret[2]+k)
end


## Strides
# lets us know if operators decouple the entries
# to split into sub problems
# A diagonal operator has essentially infinite stride
# which we represent by a factorial, so that
# the gcd with any number < 10 is the number
Base.stride(A::BandedOperator)=bandinds(A)==(0,0)?factorial(10):1
Base.stride(A::Functional)=1




## Construct operators




BandedMatrix{T<:Number}(B::Operator{T},n::Integer)=addentries!(B,bazeros(T,n,:,bandinds(B)),1:n)
BandedMatrix{T<:Number}(B::Operator{T},rws::UnitRange,::Colon)=first(rws)==1?BandedMatrix(B,last(rws)):addentries!(B,isbazeros(T,rws,:,bandinds(B)),rws).matrix


function BandedMatrix(B::Operator,kr::Range,jr::Range)
    br=bandrange(B)
    shft=kr[1]-jr[1]
    
    BandedMatrix(BandedMatrix(B,kr,:).data,length(jr),-br[1]-shft,br[end]+shft)
end


# Returns all columns in rows kr
# The first column of the returned BandedMatrix
# will be the first non-zero column

BandedMatrix(B::Operator,kr::Colon,jr::UnitRange)=BandedMatrix(B,max(1,jr[1]-bandinds(B,2)):jr[end]-bandinds(B,1),jr)



## geteindex

Base.getindex(op::Operator,k::Integer,j::Integer)=op[k:k,j:j][1,1]
Base.getindex(op::Operator,k::Integer,j::Range)=op[k:k,j][1,:]
Base.getindex(op::Operator,k::Range,j::Integer)=op[k,j:j][:,1]
Base.getindex(op::Functional,k::Integer)=op[k:k][1]

function Base.getindex(op::Functional,j::Range,k::Range)
  @assert j[1]==1 && j[end]==1
  op[k].'
end
function Base.getindex(op::Functional,j::Integer,k::Range)
  @assert j==1
  op[k].'
end



# we use slice instead of get index because we can't override
# getindex (::Colon)
# This violates the behaviour of slices though...
Base.slice(B::BandedOperator,k,j)=BandedMatrix(B,k,j)
Base.getindex(B::BandedOperator,k::Range,j::Range)=slice(B,k,j)

function subview(B::BandedOperator,kr::Range,::Colon)
     br=bandinds(B)
     BM=slice(B,kr,:)
     
     # This shifts to the correct slice
     IndexStride(BM,1-kr[1],-max(0,kr[1]-1+br[1]))
end


function subview(B::BandedOperator,::Colon,jr::Range)
     br=bandinds(B)
     BM=slice(B,:,jr)
     
     # This shifts to the correct slice
     IndexStride(BM,-max(jr[1]-1-br[end],0),1-jr[1])
end

function subview(B::BandedOperator,kr::Range,jr::Range)
     br=bandinds(B)
     BM=slice(B,kr,jr)
     
     # This shifts to the correct slice
     IndexStride(BM,1-kr[1],1-jr[1])
end



## Default addentries!
# this allows for just overriding getdiagonalentry

getdiagonalentry(B::BandedOperator,k,j)=error("Override getdiagonalentry for "*string(typeof(B)))
# 
function addentries!(B::BandedOperator,A,kr)
     br=bandinds(B)
     for k=(max(kr[1],1)):(kr[end])
         for j=max(br[1],1-k):br[end]
             A[k,k+j]=getdiagonalentry(B,k,j)
         end
     end
         
     A
end


## Default Composition with a Fun, LowRankFun, and TensorFun

Base.getindex(B::BandedOperator,f::Fun) = B*Multiplication(domainspace(B),f)
Base.getindex(B::BandedOperator,f::LowRankFun) = PlusOperator(BandedOperator[f.A[i]*B[f.B[i]] for i=1:rank(f)])
Base.getindex(B::BandedOperator,f::TensorFun) = B[LowRankFun(f)]

## Standard Operators and linear algebra


#include("ShiftOperator.jl")
include("linsolve.jl")

include("spacepromotion.jl")
include("ToeplitzOperator.jl")
include("ConstantOperator.jl")
include("TridiagonalOperator.jl")


## Operators overrided for spaces

include("Conversion.jl")
include("Multiplication.jl")
include("calculus.jl")
include("Sigma.jl")
include("Evaluation.jl")



include("SavedOperator.jl")
include("AlmostBandedOperator.jl")
include("adaptiveqr.jl")


include("algebra.jl")

include("TransposeOperator.jl")
include("StrideOperator.jl")
include("SliceOperator.jl")
include("CompactOperator.jl")


include("null.jl")
include("systems.jl")



## Conversion


Base.zero{T<:Number}(::Type{Functional{T}})=ZeroFunctional(T)
Base.zero{T<:Number}(::Type{Operator{T}})=ZeroOperator(T)
Base.zero{O<:Functional}(::Type{O})=ZeroFunctional()
Base.zero{O<:Operator}(::Type{O})=ZeroOperator()


# TODO: can convert return different type?
Base.convert{T<:Operator}(A::Type{T},n::Number)=n==0?zero(T):ConstantOperator(n)
Base.convert{T<:Operator}(A::Type{T},n::UniformScaling)=n.λ==0?zero(T):ConstantOperator(n)
Base.convert{T<:BandedOperator}(A::Type{T},f::Fun)=norm(f.coefficients)==0?zero(T):Multiplication(f)


## Promotion

for T in (:Float64,:Int64,:(Complex{Float64})), OP in (:BandedOperator,:Operator)
    @eval begin
        Base.promote_rule{N<:Number,O<:$OP{$T}}(::Type{N},::Type{O})=$OP{promote_type(N,$T)}
        Base.promote_rule{N<:Number,O<:$OP{$T}}(::Type{UniformScaling{N}},::Type{O})=$OP{promote_type(N,$T)}    
        Base.promote_rule{S,N<:Number,O<:$OP{$T}}(::Type{Fun{S,N}},::Type{O})=$OP{promote_type(N,$T)}            
    end
end

Base.promote_rule{N<:Number,O<:Operator}(::Type{N},::Type{O})=Operator
Base.promote_rule{N<:UniformScaling,O<:Operator}(::Type{N},::Type{O})=Operator
Base.promote_rule{N<:Fun,O<:Operator}(::Type{N},::Type{O})=Operator



