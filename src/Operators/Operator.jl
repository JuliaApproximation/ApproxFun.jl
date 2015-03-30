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
Base.eltype{T}(::Type{Operator{T}})=T
Base.eltype{OT<:Operator}(::Type{OT})=eltype(super(OT))




## We assume operators are T->T
rangespace(A::Operator)=AnySpace()
domainspace(A::Operator)=AnySpace()
rangespace(A::Functional)=ConstantSpace()
domain(A::Operator)=domain(domainspace(A))




Base.size(::InfiniteOperator)=[Inf,Inf]
Base.size(::Functional)=Any[1,Inf] #use Any vector so the 1 doesn't become a float
Base.size(op::Operator,k::Integer)=size(op)[k]
datalength(F::Functional)=error("Override datalength for "*string(typeof(F)))        # use datalength to indicate a finite length functional






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



bazeros{T<:Number}(B::Operator{T},n::Integer,m::Integer)=bazeros(T,n,m,bandinds(B))
bazeros{T<:Number}(B::Operator{T},n::Integer,m::Colon)=bazeros(T,n,m,bandinds(B))
bazeros{T<:Number}(B::Operator{T},n::Integer,br::(Int,Int))=bazeros(T,n,br)


BandedMatrix(B::Operator,n::Integer)=addentries!(B,bazeros(B,n,:),1:n)
BandedMatrix{T}(B::Operator{T},rws::UnitRange,::Colon)=first(rws)==1?BandedMatrix(B,last(rws)):addentries!(B,isbazeros(T,rws,:,bandinds(B)),rws).matrix

function BandedMatrix{T<:Number}(B::Operator{T},kr::StepRange,::Colon)
    stp=step(kr)

    if stp==1
        BandedMatrix(B,first(kr):last(kr),:)
    else
        str=stride(B)
        @assert mod(str,stp)==0
        # we need the shifting by bandinds to preserve mod
        @assert mod(bandinds(B,1),stp)==mod(bandinds(B,2),stp)==0
        # find column range
        jr=max(stp-mod(kr[1],stp),kr[1]+bandinds(B,1)):stp:kr[end]+bandinds(B,2)
        shf=div(first(kr)-first(jr),stp)
        bi=div(bandinds(B,1),stp)+shf,div(bandinds(B,2),stp)+shf
        A=bazeros(T,length(kr),length(jr),bi)
        addentries!(B,IndexSlice(A,first(kr)-stp,first(jr)-stp,stp,stp),kr)
        A
    end
end

function BandedMatrix(B::Operator,kr::Range,jr::Range)
    br=bandrange(B)
    shft=kr[1]-jr[1]

    BandedMatrix(BandedMatrix(B,kr,:).data,length(jr),-br[1]-shft,br[end]+shft)
end


# Returns all columns in rows kr
# The first column of the returned BandedMatrix
# will be the first non-zero column

BandedMatrix(B::Operator,kr::Colon,jr::UnitRange)=BandedMatrix(B,max(1,jr[1]-bandinds(B,2)):jr[end]-bandinds(B,1),jr)

Base.sparse(B::Operator,n::Integer)=sparse(BandedMatrix(B,n))
Base.sparse(B::Operator,n::Range,m::Range)=sparse(BandedMatrix(B,n,m))
Base.sparse(B::Operator,n::Colon,m::Range)=sparse(BandedMatrix(B,n,m))
Base.sparse(B::Operator,n::Range,m::Colon)=sparse(BandedMatrix(B,n,m))

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


## Composition with a Fun, LowRankFun, and ProductFun

Base.getindex(B::Operator,f::Fun) = B*Multiplication(domainspace(B),f)
Base.getindex{BT,S,M,T,V}(B::Operator{BT},f::LowRankFun{S,M,T,V}) = PlusOperator(BandedOperator{promote_type(BT,T,V)}[f.A[i]*B[f.B[i]] for i=1:rank(f)])
Base.getindex{BT,S,V,SS,T}(B::Operator{BT},f::ProductFun{S,V,SS,T}) = PlusOperator(BandedOperator{promote_type(BT,T)}[f.coefficients[i]*B[Fun([zeros(promote_type(BT,T),i-1),one(promote_type(BT,T))],f.space.spaces[2])] for i=1:length(f.coefficients)])

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
include("DefiniteIntegral.jl")
include("Evaluation.jl")



include("SavedOperator.jl")
include("AlmostBandedOperator.jl")
include("adaptiveqr.jl")


include("algebra.jl")

include("TransposeOperator.jl")
include("StrideOperator.jl")
include("CompactOperator.jl")


include("null.jl")
include("systems.jl")



## Conversion


Base.zero{T<:Number}(::Type{Functional{T}})=ZeroFunctional(T)
Base.zero{T<:Number}(::Type{Operator{T}})=ZeroOperator(T)
Base.zero{O<:Functional}(::Type{O})=ZeroFunctional()
Base.zero{O<:Operator}(::Type{O})=ZeroOperator()


# TODO: can convert return different type?

for TYP in (:Operator,:BandedOperator)
  @eval begin
    Base.convert{T}(A::Type{$TYP{T}},n::Number)=n==0?zero(A):ConstantOperator{T}(n)
    Base.convert{T}(A::Type{$TYP{T}},n::UniformScaling)=n.λ==0?zero(A):ConstantOperator{T}(n)
    Base.convert{T}(A::Type{$TYP{T}},f::Fun)=norm(f.coefficients)==0?zero(A):convert(A,Multiplication(f))
  end
end



## Promotion

for OP in (:BandedOperator,:Operator)
  @eval begin
      Base.promote_rule{N<:Number,O<:$OP}(::Type{N},::Type{O})=$OP{promote_type(N,eltype(O))}
      Base.promote_rule{N<:Number,O<:$OP}(::Type{UniformScaling{N}},::Type{O})=$OP{promote_type(N,eltype(O))}
      Base.promote_rule{S,N<:Number,O<:$OP}(::Type{Fun{S,N}},::Type{O})=$OP{promote_type(N,eltype(O))}
  end
end

for OP in (:BandedOperator,:Functional,:Operator)
  @eval Base.promote_rule{BO1<:$OP,BO2<:$OP}(::Type{BO1},::Type{BO2})=$OP{promote_type(eltype(BO1),eltype(BO2))}
end

