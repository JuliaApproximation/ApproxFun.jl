export Operator,Functional,InfiniteOperator
export bandinds, bandrange, linsolve, periodic
export dirichlet, neumann
export ldirichlet,rdirichlet,lneumann,rneumann
export ldiffbc,rdiffbc,diffbcs
export domainspace,rangespace


abstract Operator{T} #T is the entry type, Float64 or Complex{Float64}
abstract Functional{T} <: Operator{T}
abstract InfiniteOperator{T} <: Operator{T}   #Infinite Operators have + range
abstract BandedBelowOperator{T} <: InfiniteOperator{T}
abstract AlmostBandedOperator{T} <: BandedBelowOperator{T}
abstract BandedOperator{T} <: AlmostBandedOperator{T}

Base.eltype{T}(::Operator{T})=T
Base.eltype{T}(::Type{Operator{T}})=T
Base.eltype{OT<:Operator}(::Type{OT})=eltype(super(OT))

 #Operators are immutable
Base.copy(A::Operator)=A


## We assume operators are T->T
rangespace(A::Operator)=AnySpace()
domainspace(A::Operator)=AnySpace()
rangespace(A::Functional)=ConstantSpace()
domain(A::Operator)=domain(domainspace(A))




Base.size(::InfiniteOperator)=[Inf,Inf]
Base.size(::Functional)=Any[1,Inf] #use Any vector so the 1 doesn't become a float
Base.size(op::Operator,k::Integer)=size(op)[k]

Base.ndims(::Operator)=2
datalength(F::Functional)=Inf        # use datalength to indicate a finite length functional






## bandrange and indexrange


bandinds(A,k::Integer)=bandinds(A)[k]
bandrange(b::BandedBelowOperator)=UnitRange(bandinds(b)...)
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



bazeros(B::Operator,n::Integer,m::Integer)=bazeros(eltype(B),n,m,bandinds(B))
bazeros(B::Operator,n::Integer,m::Colon)=bazeros(eltype(B),n,m,bandinds(B))
bazeros(B::Operator,n::Integer,m::Colon,br::Tuple{Int,Int})=bazeros(eltype(B),n,m,br)


BandedMatrix(B::Operator,n::Integer)=addentries!(B,bazeros(B,n,:),1:n,:)
function BandedMatrix(B::Operator,rws::UnitRange,::Colon)
    if first(rws)==1
        BandedMatrix(B,last(rws))
    else
        addentries!(B,isbazeros(eltype(B),rws,:,bandinds(B)),rws,:).matrix
    end
end

function BandedMatrix(B::Operator,kr::StepRange,::Colon)
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
        A=bazeros(eltype(B),length(kr),length(jr),bi)
        addentries!(B,IndexSlice(A,first(kr)-stp,first(jr)-stp,stp,stp),kr,:)
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

Base.getindex(op::Operator,k::Integer,j::Range)=op[k:k,j][1,:]
Base.getindex(op::Operator,k::Range,j::Integer)=op[k,j:j][:,1]
Base.getindex(op::Functional,k::Integer)=op[k:k][1]

Base.getindex(L::BandedOperator,kr::Range,::Colon)=Functional{eltype(L)}[L[k,:] for k=kr]

function Base.getindex(op::Functional,j::Range,k::Range)
  @assert j[1]==1 && j[end]==1
  op[k].'
end
function Base.getindex(op::Functional,j::Integer,k::Range)
  @assert j==1
  op[k].'
end



## override getindex or addentries!.  Each defaults

defaultgetindex(op::Operator,k::Integer,j::Integer)=op[k:k,j:j][1,1]
defaultgetindex(B::BandedOperator,k::Range,j::Range)=slice(B,k,j)

# the defualt is to use getindex

function defaultgetindex(op::Operator,kr::Range,jr::Range)
    ret=Array(eltype(op),length(kr),length(jr))
    kk,jj=1,1
    for j=jr
        for k=kr
            ret[kk,jj]=op[k,j]
            kk+=1
        end
        kk=1
        jj+=1
    end
    ret
end


defaultgetindex(A::BandedOperator,k::Integer,::Colon)=FiniteFunctional(vec(A[k,1:1+bandinds(A,2)]),domainspace(A))
Base.getindex(B::Operator,k,j)=defaultgetindex(B,k,j)



# we use slice instead of get index because we can't override
# getindex (::Colon)
# This violates the behaviour of slices though...
Base.slice(B::BandedOperator,k,j)=BandedMatrix(B,k,j)
# Float-valued ranges are implemented to support 1:Inf to take a slice
# TODO: right now non-integer steps are only supported when the operator itself
# has compatible stride
function Base.slice(B::BandedOperator,kr::FloatRange,jr::FloatRange)
    st=step(kr)
    @assert step(jr)==st
    @assert last(kr)==last(jr)==Inf
    SliceOperator(B,first(kr)-st,first(jr)-st,st,st)
end

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
#  override either addentries! or getindex, otherwise there will be
#  an infinite loop.


function defaultaddentries!(B::BandedOperator,A,kr,::Colon)
     br=bandinds(B)
     for k=(max(kr[1],1)):(kr[end])
         for j=max(br[1],1-k):br[end]
             A[k,k+j]=B[k,k+j]
         end
     end

     A
end

addentries!(B,A,kr,::Colon)=defaultaddentries!(B,A,kr,:)


## Composition with a Fun, LowRankFun, and ProductFun

Base.getindex{BT,S,T}(B::Operator{BT},f::Fun{S,T}) = B*Multiplication(domainspace(B),f)
Base.getindex{BT,S,M,SS,T}(B::Operator{BT},f::LowRankFun{S,M,SS,T}) = mapreduce(i->f.A[i]*B[f.B[i]],+,1:rank(f))
Base.getindex{BT,S,V,SS,T}(B::Operator{BT},f::ProductFun{S,V,SS,T}) = mapreduce(i->f.coefficients[i]*B[Fun([zeros(promote_type(BT,T),i-1),one(promote_type(BT,T))],f.space[2])],+,1:length(f.coefficients))

## Standard Operators and linear algebra


include("linsolve.jl")

include("spacepromotion.jl")
include("banded/banded.jl")
include("functionals/functionals.jl")
include("almostbanded/almostbanded.jl")

include("systems.jl")

include("adaptiveqr.jl")
include("null.jl")




## Conversion


Base.zero{T<:Number}(::Type{Functional{T}})=ZeroFunctional(T)
Base.zero{T<:Number}(::Type{Operator{T}})=ZeroOperator(T)
Base.zero{O<:Functional}(::Type{O})=ZeroFunctional(eltype(O))
Base.zero{O<:Operator}(::Type{O})=ZeroOperator(eltype(O))


Base.eye(S::Space)=SpaceOperator(ConstantOperator(1.0),S,S)
Base.eye(S::Domain)=eye(Space(S))


# TODO: can convert return different type?


Base.convert{T<:Functional}(::Type{T},f::Fun)=DefiniteIntegral()[f]



Base.convert{T<:Operator}(A::Type{T},n::Number)=n==0?zero(A):ConstantOperator(eltype(T),n)
Base.convert{T<:Operator}(A::Type{T},n::UniformScaling)=n.λ==0?zero(A):ConstantOperator(eltype(T),n)
Base.convert{T<:Operator}(A::Type{T},f::Fun)=norm(f.coefficients)==0?zero(A):convert(A,Multiplication(f))





## Promotion


mat_promote_type(A,B)=promote_type(A,B)
mat_promote_type{T,B,n}(::Type{Array{T,n}},::Type{Array{B,n}})=Array{promote_type(T,B),n}
mat_promote_type{T,B<:Number,n}(::Type{Array{T,n}},::Type{B})=Array{promote_type(T,B),n}
mat_promote_type{T,B<:Number,n}(::Type{B},::Type{Array{T,n}})=Array{promote_type(T,B),n}

mat_promote_type{T,B}(::Type{BandedMatrix{T}},::Type{BandedMatrix{B}})=BandedMatrix{promote_type(T,B)}
mat_promote_type{T,B<:Number}(::Type{BandedMatrix{T}},::Type{B})=BandedMatrix{promote_type(T,B)}
mat_promote_type{T,B<:Number}(::Type{B},::Type{BandedMatrix{T}})=BandedMatrix{promote_type(T,B)}




for OP in (:BandedOperator,:Operator)
  @eval begin
      Base.promote_rule{N<:Number}(::Type{N},::Type{$OP})=$OP{N}
      Base.promote_rule{N<:Number}(::Type{UniformScaling{N}},::Type{$OP})=$OP{N}
      Base.promote_rule{S,N<:Number}(::Type{Fun{S,N}},::Type{$OP})=$OP{N}
      Base.promote_rule{N<:Number,O<:$OP}(::Type{N},::Type{O})=$OP{mat_promote_type(N,eltype(O))}
      Base.promote_rule{N<:Number,O<:$OP}(::Type{UniformScaling{N}},::Type{O})=$OP{mat_promote_type(N,eltype(O))}
      Base.promote_rule{S,N<:Number,O<:$OP}(::Type{Fun{S,N}},::Type{O})=$OP{mat_promote_type(N,eltype(O))}
  end
end

for OP in (:BandedOperator,:Functional,:Operator)
  @eval Base.promote_rule{BO1<:$OP,BO2<:$OP}(::Type{BO1},::Type{BO2})=$OP{mat_promote_type(eltype(BO1),eltype(BO2))}
end



## Wrapper

#TODO: Should cases that modify be included?
typealias WrapperOperator Union{SpaceOperator,MultiplicationWrapper,DerivativeWrapper,IntegralWrapper,
                                    ConversionWrapper,ConstantTimesOperator,ConstantTimesFunctional,TransposeOperator}
