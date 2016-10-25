export Operator
export bandinds, bandrange, linsolve, periodic
export dirichlet, neumann
export ldirichlet,rdirichlet,lneumann,rneumann
export ldiffbc,rdiffbc,diffbcs
export domainspace,rangespace


abstract Operator{T} #T is the entry type, Float64 or Complex{Float64}

Base.eltype{T}(::Operator{T}) = T
Base.eltype{T}(::Type{Operator{T}}) = T
Base.eltype{OT<:Operator}(::Type{OT}) = eltype(supertype(OT))


# default entry type
# we assume entries depend on both the domain and the basis
# realdomain case doesn't use


op_eltype(sp::Space) = promote_type(eltype(sp),prectype(domain(sp)))
op_eltype_realdomain(sp::Space) = promote_type(eltype(sp),real(prectype(domain(sp))))

 #Operators are immutable
Base.copy(A::Operator) = A


## We assume operators are T->T
rangespace(A::Operator) = error("Override rangespace for $(typeof(A))")
domainspace(A::Operator) = error("Override domainspace for $(typeof(A))")
domain(A::Operator) = domain(domainspace(A))


isconstspace(::) = false
## Functionals
isafunctional(A::Operator) = size(A,1)==1 && isconstspace(rangespace(A))


isboolvec(A) = isa(A,Repeated{Bool}) || isa(A,AbstractVector{Bool})
hastrivialblocks(A::Space) = isboolvec(blocklengths(A))
hastrivialblocks(A::Operator) = hastrivialblocks(domainspace(A)) &&
                                hastrivialblocks(rangespace(A))

macro functional(FF)
    quote
        Base.size(A::$FF,k::Integer) = k==1?1:∞
        ApproxFun.rangespace(::$FF) = ConstantSpace()
        ApproxFun.isafunctional(::$FF) = true
        ApproxFun.blockbandinds(A::$FF) = 0,hastrivialblocks(domainspace(A))?bandinds(A,2):∞
        function ApproxFun.defaultgetindex(f::$FF,k::Integer,j::Integer)
            @assert k==1
            f[j]
        end
    end
end


Base.size(A::Operator) = (size(A,1),size(A,2))
Base.size(A::Operator,k::Integer) = k==1?dimension(rangespace(A)):dimension(domainspace(A))
Base.length(A::Operator) = size(A,1) * size(A,2)


# used to compute "end" for last index
function Base.trailingsize(A::Operator, n::Integer)
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






## bandrange and indexrange
isbandedbelow(A::Operator) = isfinite(bandinds(A,1))
isbandedabove(A::Operator) = isfinite(bandinds(A,2))
isbanded(A::Operator) = isbandedbelow(A) && isbandedabove(A)


isbandedblockbandedbelow(::) = false
isbandedblockbandedabove(::) = false

isbandedblockbanded(A::Operator) = isbandedblockbandedabove(A) && isbandedblockbandedbelow(A)


# this should be determinable at compile time

#TODO: this is bad: we shouldn't assume block size 1
blockbandinds(A::Operator) = hastrivialblocks(A) && isbanded(A) ? bandinds(A) : (-∞,∞)
blockbandwidths(S::Operator) = -blockbandinds(S,1),blockbandinds(S,2)
blockbandinds(K::Operator,k::Integer) = blockbandinds(K)[k]
blockbandwidth(K::Operator,k::Integer) = k==1?-blockbandinds(K,k):blockbandinds(K,k)
subblockbandinds(K::Operator) = subblockbandinds(K,1),subblockbandinds(K,2)
subblockbandwidth(K::Operator,k::Integer) = k==1?-subblockbandinds(K,k):subblockbandinds(K,k)

isbandedblockbelow(A) = isfinite(blockbandinds(A,1))
isbandedblockabove(A) = isfinite(blockbandinds(A,2))
isbandedblock(A::Operator) = isbandedblockbelow(A) && isbandedblockabove(A)

israggedbelow(A::Operator) = isbandedbelow(A) || isbandedblockbanded(A) || isbandedblockbelow(A)


bandwidth(A::Operator) = bandwidth(A,1) + bandwidth(A,2) + 1
bandwidth(A::Operator,k::Integer) = k==1?-bandinds(A,1):bandinds(A,2)
bandwidths(A::Operator) = (bandwidth(A,1),bandwidth(A,2))
# we are always banded by the size
bandinds(A::Operator) = (1-size(A,1),size(A,2)-1)
bandinds(A,k::Integer) = bandinds(A)[k]
bandrange(b::Operator) = UnitRange(bandinds(b)...)



## Strides
# lets us know if operators decouple the entries
# to split into sub problems
# A diagonal operator has essentially infinite stride
# which we represent by a factorial, so that
# the gcd with any number < 10 is the number
Base.stride(A::Operator) =
    isdiag(A)?factorial(10):1

Base.isdiag(A::Operator) = bandinds(A)==(0,0)


## Construct operators


include("SubOperator.jl")


#
# Base.sparse(B::Operator,n::Integer)=sparse(BandedMatrix(B,n))
# Base.sparse(B::Operator,n::Range,m::Range)=sparse(BandedMatrix(B,n,m))
# Base.sparse(B::Operator,n::Colon,m::Range)=sparse(BandedMatrix(B,n,m))
# Base.sparse(B::Operator,n::Range,m::Colon)=sparse(BandedMatrix(B,n,m))

## geteindex



getindex(B::Operator,k,j) = defaultgetindex(B,k,j)
getindex(B::Operator,k) = defaultgetindex(B,k)




## override getindex.

defaultgetindex(B::Operator,k::Integer) = error("Override [k] for $(typeof(B))")
defaultgetindex(B::Operator,k::Integer,j::Integer) = error("Override [k,j] for $(typeof(B))")


# Ranges


defaultgetindex(op::Operator,kr::Range) = eltype(op)[op[k] for k in kr]
defaultgetindex(B::Operator,k::Range,j::Range) = AbstractMatrix(view(B,k,j))

defaultgetindex(op::Operator,k::Integer,j::Range) = reshape(eltype(op)[op[k,j] for j in j],1,length(j))
defaultgetindex(op::Operator,k::Range,j::Integer) = eltype(op)[op[k,j] for k in k]


defaultgetindex(A::Operator,K::Block,J::Block) = A[blockrows(A,K),blockcols(A,J)]
defaultgetindex(A::Operator,K::Block,j) = A[blockrows(A,K),j]
defaultgetindex(A::Operator,k,J::Block) = A[k,blockcols(A,J)]






# Colon casdes
defaultgetindex(A::Operator,k,j) = view(A,k,j)



# TODO: finite dimensional blocks
blockcolstop(A::Operator,K::Integer) = K-blockbandinds(A,1)
blockrowstop(A::Operator,J::Integer) = J+blockbandinds(A,2)

blockrows(A::Operator,K::Integer) = blockrange(rangespace(A),K)
blockcols(A::Operator,J::Integer) = blockrange(domainspace(A),J)


# default is to use bandwidth
# override for other shaped operators
#TODO: Why size(A,2) in colstart?
banded_colstart(A::Operator, i::Integer) = min(max(i-bandwidth(A,2), 1), size(A, 2))
banded_colstop(A::Operator, i::Integer) = min(i+bandwidth(A,1), size(A, 1))
banded_rowstart(A::Operator, i::Integer) = min(max(i-bandwidth(A,1), 1), size(A, 1))
banded_rowstop(A::Operator, i::Integer) = min(i+bandwidth(A,2), size(A, 2))

bandedblock_colstart(A::Operator, i::Integer) =
        blockstart(rangespace(A), block(domainspace(A),i)-blockbandwidth(A,2))
bandedblock_colstop(A::Operator, i::Integer) =
    min(blockstop(rangespace(A), block(domainspace(A),i)+blockbandwidth(A,1)),
        size(A, 1))
bandedblock_rowstart(A::Operator, i::Integer) =
        blockstart(domainspace(A), block(rangespace(A),i)-blockbandwidth(A,1))
bandedblock_rowstop(A::Operator, i::Integer) =
    min(blockstop(domainspace(A), block(rangespace(A),i)+blockbandwidth(A,2)),
        size(A, 2))


function bandedblockbanded_colstart(A::Operator, i::Integer)
    ds = domainspace(A)
    B = block(ds,i)
    ξ = i - blockstart(ds,B) + 1  # col in block
    bs = blockstart(rangespace(A), B-blockbandwidth(A,2))
    max(bs,bs + ξ - 1 - subblockbandwidth(A,2))
end

function bandedblockbanded_colstop(A::Operator, i::Integer)
    ds = domainspace(A)
    rs = rangespace(A)
    B = block(ds,i)
    ξ = i - blockstart(ds,B) + 1  # col in block
    Bend = B+blockbandwidth(A,1)
    bs = blockstart(rs, Bend)
    min(blockstop(rs,Bend),bs + ξ - 1 + subblockbandwidth(A,1))
end

function bandedblockbanded_rowstart(A::Operator, i::Integer)
    rs = rangespace(A)
    B = block(rs,i)
    ξ = i - blockstart(rs,B) + 1  # row in block
    bs = blockstart(domainspace(A), B-blockbandwidth(A,1))
    max(bs,bs + ξ - 1 - subblockbandwidth(A,1))
end

function bandedblockbanded_rowstop(A::Operator, i::Integer)
    ds = domainspace(A)
    rs = rangespace(A)
    B = block(rs,i)
    ξ = i - blockstart(rs,B) + 1  # row in block
    Bend = B+blockbandwidth(A,2)
    bs = blockstart(ds, Bend)
    min(blockstop(ds,Bend),bs + ξ - 1 + subblockbandwidth(A,2))
end


unstructured_colstart(A, i) = 1
unstructured_colstop(A, i) = size(A,1)
unstructured_rowstart(A, i) = 1
unstructured_rowstop(A, i) = size(A,2)


function default_colstart(A::Operator, i::Integer)
    if isbandedabove(A)
        banded_colstart(A,i)
    elseif isbandedblockbanded(A)
        bandedblockbanded_colstart(A, i)
    elseif isbandedblock(A)
        bandedblock_colstart(A, i)
    else
        unstructured_colstart(A, i)
    end
end

function default_colstop(A::Operator, i::Integer)
    if isbandedbelow(A)
        banded_colstop(A,i)
    elseif isbandedblockbanded(A)
        bandedblockbanded_colstop(A, i)
    elseif isbandedblock(A)
        bandedblock_colstop(A, i)
    else
        unstructured_colstop(A, i)
    end
end

function default_rowstart(A::Operator, i::Integer)
    if isbandedbelow(A)
        banded_rowstart(A,i)
    elseif isbandedblockbanded(A)
        bandedblockbanded_rowstart(A, i)
    elseif isbandedblock(A)
        bandedblock_rowstart(A, i)
    else
        unstructured_rowstart(A, i)
    end
end

function default_rowstop(A::Operator, i::Integer)
    if isbandedabove(A)
        banded_rowstop(A,i)
    elseif isbandedblockbanded(A)
        bandedblockbanded_rowstop(A, i)
    elseif isbandedblock(A)
        bandedblock_rowstop(A, i)
    else
        unstructured_rowstop(A, i)
    end
end



for OP in (:colstart,:colstop,:rowstart,:rowstop)
    defOP = parse("default_"*string(OP))
    @eval begin
        $OP(A::Operator,i::Integer) = $defOP(A,i)
        $OP(A::Operator,i::Infinity{Bool}) = ∞
    end
end




function defaultgetindex(A::Operator,::Type{FiniteRange},::Type{FiniteRange})
    if isfinite(size(A,1)) && isfinite(size(A,2))
        A[1:size(A,1),1:size(A,2)]
    else
        error("Only exists for finite operators.")
    end
end

defaultgetindex(A::Operator,k::Type{FiniteRange},J::Block) = A[k,blockcols(A,J)]
function defaultgetindex(A::Operator,::Type{FiniteRange},jr)
    cs=mapreduce(j->colstop(A,j),max,jr)
    A[1:cs,jr]
end

defaultgetindex(A::Operator,K::Block,j::Type{FiniteRange}) = A[blockrows(A,K),j]
defaultgetindex(A::Operator,kr,::Type{FiniteRange}) =
    A[kr,1:rowstop(A,maximum(kr))]





## Composition with a Fun, LowRankFun, and ProductFun

defaultgetindex{BT,S,T}(B::Operator{BT},f::Fun{S,T}) = B*Multiplication(domainspace(B),f)
defaultgetindex{BT,S,M,SS,T}(B::Operator{BT},f::LowRankFun{S,M,SS,T}) =
    mapreduce(i->f.A[i]*B[f.B[i]],+,1:rank(f))
defaultgetindex{BT,S,V,SS,T}(B::Operator{BT},f::ProductFun{S,V,SS,T}) =
    mapreduce(i->f.coefficients[i]*B[Fun([zeros(promote_type(BT,T),i-1);
                                            one(promote_type(BT,T))],f.space[2])],
                +,1:length(f.coefficients))



# Convenience for wrapper ops
unwrap_axpy!(α,P,A) = BLAS.axpy!(α,view(parent(P).op,P.indexes[1],P.indexes[2]),A)
iswrapper(::)=false
haswrapperstructure(::)=false

# use this for wrapper operators that have the same structure but
# not necessarily the same entries
#
#  Ex: c*op or real(op)
macro wrapperstructure(Wrap)
    ret = quote
        iswrapperstructure(::$Wrap) = true
    end

    for func in (:(ApproxFun.bandinds),:(Base.stride),
                 :(ApproxFun.isbandedblockbanded),:(ApproxFun.isbandedblock),
                 :(ApproxFun.israggedbelow),:(Base.size),:(ApproxFun.isbanded),
                 :(ApproxFun.bandwidth),:(ApproxFun.bandwidths))
        ret = quote
            $ret

            $func(D::$Wrap) = $func(D.op)
        end
    end

     for func in (:(ApproxFun.bandwidth),:(ApproxFun.colstart),:(ApproxFun.colstop),
                     :(ApproxFun.rowstart),:(ApproxFun.rowstop),:(ApproxFun.blockbandinds),
                     :(Base.size),:(ApproxFun.bandinds))
         ret = quote
             $ret

             $func(D::$Wrap,k::Integer) = $func(D.op,k)
         end
     end

    esc(ret)
end



# use this for wrapper operators that have the same entries but
# not necessarily the same spaces
#
macro wrappergetindex(Wrap)
    ret = quote
        Base.getindex(OP::$Wrap,k::Integer...) =
            OP.op[k...]::eltype(OP)

        BLAS.axpy!{T,OP<:$Wrap}(α,P::ApproxFun.SubOperator{T,OP},A::AbstractMatrix) =
            ApproxFun.unwrap_axpy!(α,P,A)
    end

    for TYP in (:(BandedMatrices.BandedMatrix),:(ApproxFun.BandedBlockBandedMatrix),
                :Matrix,:AbstractMatrix,:Vector,:AbstractVector)
        ret = quote
            $ret

            Base.convert{T,OP<:$Wrap}(::Type{$TYP},P::ApproxFun.SubOperator{T,OP}) =
                $TYP(view(parent(P).op,P.indexes[1],P.indexes[2]))
        end
    end

    ret = quote
        $ret

        ApproxFun.@wrapperstructure($Wrap) # structure is automatically inherited
    end

    esc(ret)
end

# use this for wrapper operators that have the same spaces but
# not necessarily the same entries or structure
#
macro wrapperspaces(Wrap)
    ret = quote  end

    for func in (:(ApproxFun.rangespace),:(ApproxFun.domain),
                 :(ApproxFun.domainspace),:(ApproxFun.isconstop))
        ret = quote
            $ret

            $func(D::$Wrap) = $func(D.op)
        end
    end

    esc(ret)
end


# use this for wrapper operators that have the same entries and same spaces
#
macro wrapper(Wrap)
    ret = quote
        ApproxFun.@wrappergetindex($Wrap)
        ApproxFun.@wrapperspaces($Wrap)

        ApproxFun.iswrapper(::$Wrap) = true
    end


    esc(ret)
end

## Standard Operators and linear algebra



include("linsolve.jl")

include("spacepromotion.jl")
include("banded/banded.jl")
include("general/general.jl")

include("functionals/functionals.jl")
include("almostbanded/almostbanded.jl")

include("systems.jl")

include("qrfact.jl")
include("nullspace.jl")




## Conversion



Base.zero{T<:Number}(::Type{Operator{T}}) = ZeroOperator(T)
Base.zero{O<:Operator}(::Type{O}) = ZeroOperator(eltype(O))


Base.eye(S::Space) = IdentityOperator(S)
Base.eye(S::Domain) = eye(Space(S))

Base.convert{T}(A::Type{Operator{T}},f::Fun) =
    norm(f.coefficients)==0?zero(A):convert(A,Multiplication(f))





## Promotion





Base.promote_rule{N<:Number}(::Type{N},::Type{Operator}) = Operator{N}
Base.promote_rule{N<:Number}(::Type{UniformScaling{N}},::Type{Operator}) =
    Operator{N}
Base.promote_rule{S,N<:Number}(::Type{Fun{S,N}},::Type{Operator}) = Operator{N}
Base.promote_rule{N<:Number,O<:Operator}(::Type{N},::Type{O}) =
    Operator{promote_type(N,eltype(O))}
Base.promote_rule{N<:Number,O<:Operator}(::Type{UniformScaling{N}},::Type{O}) =
    Operator{promote_type(N,eltype(O))}
Base.promote_rule{S,N<:Number,O<:Operator}(::Type{Fun{S,N}},::Type{O}) =
    Operator{promote_type(N,eltype(O))}

Base.promote_rule{BO1<:Operator,BO2<:Operator}(::Type{BO1},::Type{BO2}) =
    Operator{promote_type(eltype(BO1),eltype(BO2))}




## Wrapper

#TODO: Should cases that modify be included?
typealias WrapperOperator Union{SpaceOperator,MultiplicationWrapper,DerivativeWrapper,IntegralWrapper,
                                    ConversionWrapper,ConstantTimesOperator,TransposeOperator}





# The following support converting an Operator to a Matrix or BandedMatrix

## BLAS and matrix routines
# We assume that copy may be overriden

BLAS.axpy!(a,X::Operator,Y::AbstractMatrix) = BLAS.axpy!(a,AbstractMatrix(X),Y)

# this is for operators that implement copy via axpy!

bzeros(S::Operator) = bzeros(eltype(S),size(S,1),size(S,2),bandwidth(S,1),bandwidth(S,2))
Base.zeros(S::Operator) = zeros(eltype(S),size(S,1),size(S,2))
bbbzeros(S::Operator) = bbbzeros(eltype(S),blockbandwidth(S,1),blockbandwidth(S,2),
                    subblockbandwidth(S,1),subblockbandwidth(S,2),
            blocklengthrange(rangetensorizer(S),1:size(S,1)),
            blocklengthrange(domaintensorizer(S),1:size(S,2)))

rzeros(S::Operator) = rzeros(eltype(S),size(S,1),Int[max(0,colstop(S,j)) for j=1:size(S,2)])

for (TYP,ZERS) in ((:BandedMatrix,:bzeros),(:Matrix,:zeros),
                   (:BandedBlockBandedMatrix,:bbbzeros),
                   (:RaggedMatrix,:rzeros),(:BandedBlockMatrix,:bbzeros))
    @eval convert_axpy!(::Type{$TYP},S::Operator) =
        BLAS.axpy!(one(eltype(S)),S,$ZERS(S))
end




function Base.convert(::Type{Matrix},S::Operator)
   if isinf(size(S,1)) || isinf(size(S,2))
       error("Cannot convert $S to a Matrix")
   end

   eltype(S)[S[k,j] for k=1:size(S,1),j=1:size(S,2)]
end

Base.convert(::Type{BandedMatrix},S::Operator) = default_bandedmatrix(S)
function Base.convert(::Type{RaggedMatrix},S::Operator)
    if isbanded(S)
        RaggedMatrix(BandedMatrix(S))
    elseif isbandedblockbanded(S)
        RaggedMatrix(BandedBlockBandedMatrix(S))
    else
        default_raggedmatrix(S)
    end
end

function Base.convert(::Type{Vector},S::Operator)
    if size(S,2) ≠ 1  || isinf(size(S,1))
        error("Cannot convert $S to a AbstractVector")
    end

    eltype(S)[S[k] for k=1:size(S,1)]
end



Base.convert(::Type{AbstractMatrix},S::Operator) = Matrix(S)

function Base.convert(::Type{AbstractMatrix},S::SubOperator)
    if isinf(size(S,1)) || isinf(size(S,2))
        throw(BoundsError())
    end
    if isbanded(parent(S))
        BandedMatrix(S)
    elseif isbandedblockbanded(parent(S))
        BandedBlockBandedMatrix(S)
    elseif israggedbelow(parent(S))
        RaggedMatrix(S)
    else
        Matrix(S)
    end
end

Base.convert(::Type{AbstractVector},S::Operator) = Vector(S)
