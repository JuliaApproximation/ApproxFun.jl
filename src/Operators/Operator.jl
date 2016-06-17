export Operator,BandedOperator,Functional,InfiniteOperator
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


# default entry type
# we assume entries depend on both the domain and the basis
# realdomain case doesn't use


op_eltype(sp::Space)=promote_type(eltype(sp),eltype(domain(sp)))
op_eltype_realdomain(sp::Space)=promote_type(eltype(sp),real(eltype(domain(sp))))

 #Operators are immutable
Base.copy(A::Operator)=A


## We assume operators are T->T
rangespace(A::Operator)=AnySpace()
domainspace(A::Operator)=AnySpace()
rangespace(A::Functional)=ConstantSpace()
domain(A::Operator)=domain(domainspace(A))




Base.size(A::Operator) = (size(A,1),size(A,2))
Base.size(A::Operator,k::Integer) = k==1?dimension(rangespace(A)):dimension(domainspace(A))

# used to compute "end" for last index
function Base.trailingsize(A::BandedOperator, n::Integer)
    if n > 2
        1
    elseif n==2
        size(A,2)
    elseif isinf(size(A,2)) || isinf(size(A,1))
        ∞
    else
        size(A,1)*size(A,2)
    end
end

Base.ndims(::Operator) = 2
datalength(F::Functional) = ∞        # use datalength to indicate a finite length functional






## bandrange and indexrange

bandwidth(A::BandedOperator,k::Integer)=k==1?-bandinds(A,1):bandinds(A,2)
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
Base.stride(A::BandedOperator)=isdiag(A)?factorial(10):1
Base.stride(A::Functional)=1

Base.isdiag(A::BandedOperator)=bandinds(A)==(0,0)


## Construct operators


include("SubMatrix.jl")



bzeros(B::Operator,n::Integer,m::Integer)=bzeros(eltype(B),n,m,bandinds(B))
bzeros(B::Operator,n::Integer,m::Colon)=bzeros(eltype(B),n,m,bandinds(B))
bzeros(B::Operator,n::Integer,m::Colon,br::Tuple{Int,Int})=bzeros(eltype(B),n,m,br)


# The following support converting an Operator to a BandedMatrix
#  In this case : is interpreted to mean all nonzero rows

BandedMatrix(B::Operator,kr::Range,jr::Range) =
    copy(sub(B,kr,jr))

BandedMatrix(B::Operator,rws::Range,::Colon) =
    BandedMatrix(B,rws,1:last(rws)+bandwidth(B,2))

BandedMatrix(B::Operator,kr::Colon,jr::UnitRange) =
    BandedMatrix(B,max(1,jr[1]-bandinds(B,2)):jr[end]-bandinds(B,1),jr)
#
# Base.sparse(B::Operator,n::Integer)=sparse(BandedMatrix(B,n))
# Base.sparse(B::Operator,n::Range,m::Range)=sparse(BandedMatrix(B,n,m))
# Base.sparse(B::Operator,n::Colon,m::Range)=sparse(BandedMatrix(B,n,m))
# Base.sparse(B::Operator,n::Range,m::Colon)=sparse(BandedMatrix(B,n,m))

## geteindex



Base.getindex(B::Operator,k,j) = defaultgetindex(B,k,j)
Base.getindex(B::Operator,k) = defaultgetindex(B,k)


## override getindex.

defaultgetindex(B::Operator,k::Integer) = error("Override getindex for $(typeof(B))")
defaultgetindex(B::Operator,k::Integer,j::Integer) = error("Override getindex for $(typeof(B))")


# Ranges


defaultgetindex(op::Operator,kr::Range)=eltype(op)[op[k] for k in kr]
defaultgetindex(B::Operator,k::Range,j::Range) = copy(sub(B,k,j))

defaultgetindex(op::Operator,k::Integer,j::Range) = reshape(eltype(op)[op[k,j] for j in j],1,length(j))
defaultgetindex(op::Operator,k::Range,j::Integer) = eltype(op)[op[k,j] for k in k]


function defaultgetindex(op::Functional,k::Integer,j::Integer)
    @assert k==1
    op[j]
end







# Colon casdes
defaultgetindex(L::BandedOperator,kr::Range,::Colon)=Functional{eltype(L)}[L[k,:] for k=kr]
defaultgetindex(A::BandedOperator,k::Integer,::Colon) =
    FiniteFunctional(vec(A[k,1:1+bandinds(A,2)]),domainspace(A))
defaultgetindex(A::Operator,kr::Range,::Colon) = sub(A,kr,:)
defaultgetindex(A::Operator,::Colon,jr::Range) = sub(A,:,jr)
defaultgetindex(A::Operator,::Colon,::Colon) = A

defaultgetindex(A::Operator,kr::AbstractCount,jr::AbstractCount) = sub(A,kr,jr)
defaultgetindex(B::Operator,k::AbstractCount,::Colon) = B[k,1:end]
defaultgetindex(B::Operator,::Colon,j::AbstractCount) = B[1:end,j]




## Composition with a Fun, LowRankFun, and ProductFun

defaultgetindex{BT,S,T}(B::Operator{BT},f::Fun{S,T}) = B*Multiplication(domainspace(B),f)
defaultgetindex{BT,S,M,SS,T}(B::Operator{BT},f::LowRankFun{S,M,SS,T}) = mapreduce(i->f.A[i]*B[f.B[i]],+,1:rank(f))
defaultgetindex{BT,S,V,SS,T}(B::Operator{BT},f::ProductFun{S,V,SS,T}) =
    mapreduce(i->f.coefficients[i]*B[Fun([zeros(promote_type(BT,T),i-1);one(promote_type(BT,T))],f.space[2])],+,1:length(f.coefficients))



# Convenience for wrapper ops
unwrap_axpy!(α,P,A) = BLAS.axpy!(α,sub(parent(P).op,P.indexes[1],P.indexes[2]),A)
iswrapper(::)=false


macro wrappergetindex(Wrap)
    ret = quote
        Base.getindex(OP::$Wrap,k::Integer,j::Integer) =
            OP.op[k,j]

        BLAS.axpy!{T,OP<:$Wrap}(α,P::ApproxFun.SubBandedMatrix{T,OP},A::AbstractMatrix) =
            ApproxFun.unwrap_axpy!(α,P,A)

        Base.copy{T,OP<:$Wrap}(P::ApproxFun.SubBandedMatrix{T,OP}) =
            copy(sub(parent(P).op,P.indexes[1],P.indexes[2]))
    end

    esc(ret)
end


macro wrapper(Wrap)
    ret = quote
        @wrappergetindex($Wrap)

        ApproxFun.iswrapper(::$Wrap)=true
    end
    for func in (:(ApproxFun.rangespace),:(ApproxFun.domainspace),
                 :(ApproxFun.bandinds),:(ApproxFun.domain),:(Base.stride))
        ret=quote
            $ret

            $func(D::$Wrap)=$func(D.op)
        end
    end
    esc(ret)
end

## Standard Operators and linear algebra



include("linsolve.jl")

include("spacepromotion.jl")
include("banded/banded.jl")
include("functionals/functionals.jl")
include("almostbanded/almostbanded.jl")

include("systems.jl")

include("adaptiveqr.jl")
include("nullspace.jl")




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
