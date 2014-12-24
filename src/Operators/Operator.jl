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


#TODO: Speed up by taking only slice
Base.getindex(B::BandedOperator,k::Range,j::Range)=BandedMatrix(B,max(k[end],j[end]))[k,j]



## bandrange and indexrange


bandinds(A::BandedOperator,k::Integer)=bandinds(A)[k]
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


ShiftMatrix{T<:Number}(B::Operator{T},n::Integer)=addentries!(B,sazeros(T,n,bandinds(B)),1:n)
ShiftMatrix{T<:Number}(B::Operator{T},rws::Range)=addentries!(B,issazeros(T,rws,bandinds(B)),rws).matrix
ShiftMatrix{T<:Number}(B::Operator{T},rws::(Int,Int))=addentries!(B,issazeros(T,rws,bandinds(B)),rws[1]:rws[end]).matrix
BandedMatrix(B::Operator,n)=BandedMatrix(ShiftMatrix(B,n))
BandedMatrix(B::Operator,::Colon,col::Integer)=BandedMatrix(ShiftMatrix(B,col-bandinds(B,1)),col)


ShiftArray{T<:Number}(B::Operator{T},k::Range,j::Range)=addentries!(B,sazeros(T,k,j),k)
ShiftArray(B::Operator,k::Range)=ShiftArray(B,k,bandrange(B))
BandedArray(B::Operator,k::Range)=BandedArray(B,k,(k[1]+bandinds(B)[1]):(k[end]+bandinds(B)[end]))
BandedArray(B::Operator,k::Range,cs)=BandedArray(ShiftArray(B,k,bandrange(B)),cs)


## Default addentries!
# this allows for just overriding getdiagonalentry

getdiagonalentry(B::BandedOperator,k,j)=error("Override either getdiagonalentry or addentries! for "*string(typeof(B)))

function addentries!(B::BandedOperator,A,kr)
    if isa(A,BandedMatrix)
        A=BandedMatrix(addentries!(B,ShiftMatrix(A),kr))
    else
        br=bandinds(B)
        for k=(max(kr[1],1)):(kr[end])
            for j=max(br[1],1-k):br[end]
                A[k,j]=getdiagonalentry(B,k,j)
            end
        end
    end
        
    A
end


## Default Composition with a Fun, LowRankFun, and TensorFun

Base.getindex(B::BandedOperator,f::Fun) = B*Multiplication(domainspace(B),f)
Base.getindex(B::BandedOperator,f::LowRankFun) = PlusOperator(BandedOperator[f.A[i]*B[f.B[i]] for i=1:rank(f)])
Base.getindex(B::BandedOperator,f::TensorFun) = B[LowRankFun(f)]

## Standard Operators and linear algebra


include("ShiftOperator.jl")
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



