export Operator
export bandinds, bandrange, \, periodic
export dirichlet, neumann
export ldirichlet,rdirichlet,lneumann,rneumann
export ldiffbc,rdiffbc,diffbcs
export domainspace,rangespace


abstract type Operator{T} end #T is the entry type, Float64 or Complex{Float64}

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
# block lengths of a space are 1
hastrivialblocks(A::Space) = isboolvec(blocklengths(A))
hastrivialblocks(A::Operator) = hastrivialblocks(domainspace(A)) &&
                                hastrivialblocks(rangespace(A))

# blocklengths are constant lengths
hasconstblocks(A::Space) = isa(blocklengths(A),Repeated)
hasconstblocks(A::Operator) = hasconstblocks(domainspace(A)) && hasconstblocks(rangespace(A)) &&
                                blocklengths(domainspace(A)).x == blocklengths(rangespace(A)).x


macro functional(FF)
    quote
        Base.size(A::$FF,k::Integer) = k==1?1:∞
        ApproxFun.rangespace(::$FF) = ConstantSpace()
        ApproxFun.isafunctional(::$FF) = true
        ApproxFun.blockbandinds(A::$FF) = 0,hastrivialblocks(domainspace(A))?bandinds(A,2):∞
        function ApproxFun.defaultgetindex(f::$FF,k::Integer,j::Integer)
            @assert k==1
            f[j]::eltype(f)
        end
    end
end


blocksize(A::Operator,k) = k==1 ? length(blocklengths(rangespace(A))) : length(blocklengths(domainspace(A)))
blocksize(A::Operator) = (blocksize(A,1),blocksize(A,2))


Base.size(A::Operator) = (size(A,1),size(A,2))
Base.size(A::Operator,k::Integer) = k==1 ? dimension(rangespace(A)) : dimension(domainspace(A))
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
#TODO: I think it can be generalized to the case when the domainspace
# blocklengths == rangespace blocklengths, in which case replace the definition
# of p with maximum(blocklength(domainspace(A)))
function blockbandinds(A::Operator)
    hastrivialblocks(A) && return bandinds(A)

    if hasconstblocks(A)
        a,b = bandinds(A)
        p = blocklengths(domainspace(A)).x
        return (fld(a,p),-fld(-b,p))
    end

    #TODO: Generalize to finite dimensional
    if size(A,2) == 1
        rs = rangespace(A)

        if hasconstblocks(rs)
            a = bandinds(A,1)
            p = blocklengths(rs).x
            return (fld(a,p),0)
        end
    end


    return (-∞,∞)
end

# assume dense blocks
subblockbandinds(K::Operator,k) = k==1 ? 1-maximum(blocklengths(rangespace(K))) : maximum(blocklengths(domainspace(K)))-1

isblockbandedbelow(A) = isfinite(blockbandinds(A,1))
isblockbandedabove(A) = isfinite(blockbandinds(A,2))
isblockbanded(A::Operator) = isblockbandedbelow(A) && isblockbandedabove(A)

israggedbelow(A::Operator) = isbandedbelow(A) || isbandedblockbanded(A) || isblockbandedbelow(A)


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
defaultgetindex(B::Operator,k::Block,j::Block) = AbstractMatrix(view(B,k,j))
defaultgetindex(B::Operator,k::Range,j::Block) = AbstractMatrix(view(B,k,j))
defaultgetindex(B::Operator,k::Block,j::Range) = AbstractMatrix(view(B,k,j))
defaultgetindex(B::Operator,k::Range,j::Range) = AbstractMatrix(view(B,k,j))

defaultgetindex(op::Operator,k::Integer,j::Range) = eltype(op)[op[k,j] for j in j]
defaultgetindex(op::Operator,k::Range,j::Integer) = eltype(op)[op[k,j] for k in k]


# Colon casdes
defaultgetindex(A::Operator,kj::CartesianIndex{2}) = A[kj[1],kj[2]]
defaultgetindex(A::Operator,kj::CartesianIndex{1}) = A[kj[1]]
defaultgetindex(A::Operator,k,j) = view(A,k,j)



# TODO: finite dimensional blocks
blockcolstart(A::Operator,J::Integer) = Block(max(1,J-blockbandwidth(A,2)))
blockrowstart(A::Operator,K::Integer) = Block(max(1,K-blockbandwidth(A,1)))
blockcolstop(A::Operator,J::Integer) = Block(min(J+blockbandwidth(A,1),blocksize(A,1)))
blockrowstop(A::Operator,K::Integer) = Block(min(K+blockbandwidth(A,2),blocksize(A,2)))

blockrows(A::Operator,K::Integer) = blockrange(rangespace(A),K)
blockcols(A::Operator,J::Integer) = blockrange(domainspace(A),J)


# default is to use bandwidth
# override for other shaped operators
#TODO: Why size(A,2) in colstart?
banded_colstart(A::Operator, i::Integer) = min(max(i-bandwidth(A,2), 1), size(A, 2))
banded_colstop(A::Operator, i::Integer) = max(0,min(i+bandwidth(A,1), size(A, 1)))
banded_rowstart(A::Operator, i::Integer) = min(max(i-bandwidth(A,1), 1), size(A, 1))
banded_rowstop(A::Operator, i::Integer) = max(0,min(i+bandwidth(A,2), size(A, 2)))

blockbanded_colstart(A::Operator, i::Integer) =
        blockstart(rangespace(A), block(domainspace(A),i)-blockbandwidth(A,2))
blockbanded_colstop(A::Operator, i::Integer) =
    min(blockstop(rangespace(A), block(domainspace(A),i)+blockbandwidth(A,1)),
        size(A, 1))
blockbanded_rowstart(A::Operator, i::Integer) =
        blockstart(domainspace(A), block(rangespace(A),i)-blockbandwidth(A,1))
blockbanded_rowstop(A::Operator, i::Integer) =
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
    elseif isblockbanded(A)
        blockbanded_colstart(A, i)
    else
        unstructured_colstart(A, i)
    end
end

function default_colstop(A::Operator, i::Integer)
    if isbandedbelow(A)
        banded_colstop(A,i)
    elseif isbandedblockbanded(A)
        bandedblockbanded_colstop(A, i)
    elseif isblockbanded(A)
        blockbanded_colstop(A, i)
    else
        unstructured_colstop(A, i)
    end
end

function default_rowstart(A::Operator, i::Integer)
    if isbandedbelow(A)
        banded_rowstart(A,i)
    elseif isbandedblockbanded(A)
        bandedblockbanded_rowstart(A, i)
    elseif isblockbanded(A)
        blockbanded_rowstart(A, i)
    else
        unstructured_rowstart(A, i)
    end
end

function default_rowstop(A::Operator, i::Integer)
    if isbandedabove(A)
        banded_rowstop(A,i)
    elseif isbandedblockbanded(A)
        bandedblockbanded_rowstop(A, i)
    elseif isblockbanded(A)
        blockbanded_rowstop(A, i)
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
function defaultgetindex(A::Operator,::Type{FiniteRange},jr::AbstractVector{Int})
    cs = (isbanded(A) || isblockbandedbelow(A)) ? colstop(A,maximum(jr)) : mapreduce(j->colstop(A,j),max,jr)
    A[1:cs,jr]
end

function defaultgetindex(A::Operator,::Type{FiniteRange},jr::AbstractVector{Block})
    cs = (isbanded(A) || isblockbandedbelow(A)) ? blockcolstop(A,maximum(jr)) : mapreduce(j->blockcolstop(A,j),max,jr)
    A[Block(1):cs,jr]
end

function Base.view(A::Operator,::Type{FiniteRange},jr::AbstractVector{Int})
    cs = (isbanded(A) || isblockbandedbelow(A)) ? colstop(A,maximum(jr)) : mapreduce(j->colstop(A,j),max,jr)
    view(A,1:cs,jr)
end

function Base.view(A::Operator,::Type{FiniteRange},jr::AbstractVector{Block})
    cs = (isbanded(A) || isblockbandedbelow(A)) ? blockcolstop(A,maximum(jr)) : mapreduce(j->blockcolstop(A,j),max,jr)
    view(A,Block(1):cs,jr)
end


defaultgetindex(A::Operator,K::Block,j::Type{FiniteRange}) = A[blockrows(A,K),j]
defaultgetindex(A::Operator,kr,::Type{FiniteRange}) =
    A[kr,1:rowstop(A,maximum(kr))]





## Composition with a Fun, LowRankFun, and ProductFun

defaultgetindex(B::Operator,f::Fun) = B*Multiplication(domainspace(B),f)
defaultgetindex(B::Operator,f::LowRankFun) = mapreduce(i->f.A[i]*B[f.B[i]],+,1:rank(f))
defaultgetindex{BT,S,V,SS,T}(B::Operator{BT},f::ProductFun{S,V,SS,T}) =
    mapreduce(i->f.coefficients[i]*B[Fun(f.space[2],[zeros(promote_type(BT,T),i-1);
                                            one(promote_type(BT,T))])],
                +,1:length(f.coefficients))



# Convenience for wrapper ops
unwrap_axpy!(α,P,A) = BLAS.axpy!(α,view(parent(P).op,P.indexes[1],P.indexes[2]),A)
iswrapper(::) = false
haswrapperstructure(::) = false

# use this for wrapper operators that have the same structure but
# not necessarily the same entries
#
#  Ex: c*op or real(op)
macro wrapperstructure(Wrap)
    ret = quote
        haswrapperstructure(::$Wrap) = true
    end

    for func in (:(ApproxFun.bandinds),:(Base.stride),
                 :(ApproxFun.isbandedblockbanded),:(ApproxFun.isblockbanded),
                 :(ApproxFun.israggedbelow),:(Base.size),:(ApproxFun.isbanded),
                 :(ApproxFun.bandwidth),:(ApproxFun.bandwidths),
                 :(ApproxFun.blockbandinds),:(ApproxFun.subblockbandinds),
                 :(Base.issymmetric))
        ret = quote
            $ret

            $func(D::$Wrap) = $func(D.op)
        end
    end

     for func in (:(ApproxFun.bandwidth),:(ApproxFun.colstart),:(ApproxFun.colstop),
                     :(ApproxFun.rowstart),:(ApproxFun.rowstop),:(ApproxFun.blockbandinds),
                     :(Base.size),:(ApproxFun.bandinds),:(ApproxFun.subblockbandinds))
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

        Base.getindex(OP::$Wrap,k...) = OP.op[k...]

        BLAS.axpy!{T,OP<:$Wrap}(α,P::ApproxFun.SubOperator{T,OP},A::AbstractMatrix) =
            ApproxFun.unwrap_axpy!(α,P,A)

        A_mul_B_coefficients(A::$Wrap,b) = A_mul_B_coefficients(A.op,b)
        A_mul_B_coefficients{T,OP<:$Wrap}(A::ApproxFun.SubOperator{T,OP,Tuple{UnitRange{Int},UnitRange{Int}}},b) =
            A_mul_B_coefficients(view(parent(A).op,S.indexes[1],S.indexes[2]),b)
        A_mul_B_coefficients{T,OP<:$Wrap}(A::ApproxFun.SubOperator{T,OP},b) =
            A_mul_B_coefficients(view(parent(A).op,S.indexes[1],S.indexes[2]),b)
    end

    for TYP in (:(BandedMatrices.BandedMatrix),:Matrix,:Vector,:AbstractVector)
        ret = quote
            $ret

            Base.convert{T,OP<:$Wrap}(::Type{$TYP},P::ApproxFun.SubOperator{T,OP}) =
                $TYP(view(parent(P).op,P.indexes[1],P.indexes[2]))
        end
    end

    ret = quote
        $ret

        # fast converts to banded matrices would be based on indices, not blocks
        function Base.convert{T,OP<:$Wrap}(::Type{BandedMatrices.BandedMatrix},
                                S::ApproxFun.SubOperator{T,OP,Tuple{UnitRange{ApproxFun.Block},UnitRange{ApproxFun.Block}}})
            A = parent(S)
            ds = domainspace(A)
            rs = rangespace(A)
            KR,JR = parentindexes(S)
            BandedMatrix(view(A,
                              blockstart(rs,KR[1]):blockstop(rs,KR[end]),
                              blockstart(ds,JR[1]):blockstop(ds,JR[end])))
        end


        # if the spaces change, then we need to be smarter
        function Base.convert{T,OP<:$Wrap}(::Type{ApproxFun.BlockBandedMatrix},S::ApproxFun.SubOperator{T,OP})
            P = parent(S)
            if blocklengths(domainspace(P)) == blocklengths(domainspace(P.op)) &&
                    blocklengths(rangespace(P)) == blocklengths(rangespace(P.op))
                BlockBandedMatrix(view(parent(S).op,S.indexes[1],S.indexes[2]))
            else
                default_blockbandedmatrix(S)
            end
        end

        function Base.convert{T,OP<:$Wrap}(::Type{ApproxFun.BandedBlockBandedMatrix},S::ApproxFun.SubOperator{T,OP})
            P = parent(S)
            if blocklengths(domainspace(P)) == blocklengths(domainspace(P.op)) &&
                    blocklengths(rangespace(P)) == blocklengths(rangespace(P.op))
                BandedBlockBandedMatrix(view(parent(S).op,S.indexes[1],S.indexes[2]))
            else
                default_bandedblockbandedmatrix(S)
            end
        end

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



include("A_ldiv_B.jl")

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

Base.convert(A::Type{Operator},f::Fun) =
    norm(f.coefficients)==0?ZeroOperator():Multiplication(f)





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
const WrapperOperator = Union{SpaceOperator,MultiplicationWrapper,DerivativeWrapper,IntegralWrapper,
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
bbzeros(S::Operator) = bbzeros(eltype(S),blockbandwidth(S,1),blockbandwidth(S,2),
            blocklengths(rangespace(S)),blocklengths(domainspace(S)))

rzeros(S::Operator) = rzeros(eltype(S),size(S,1),Int[max(0,colstop(S,j)) for j=1:size(S,2)])

for (TYP,ZERS) in ((:BandedMatrix,:bzeros),(:Matrix,:zeros),
                   (:BandedBlockBandedMatrix,:bbbzeros),
                   (:RaggedMatrix,:rzeros),(:BlockBandedMatrix,:bbzeros))
    @eval convert_axpy!(::Type{$TYP},S::Operator) =
        BLAS.axpy!(one(eltype(S)),S,$ZERS(S))
end




function Base.convert(::Type{Matrix},S::Operator)
   if isinf(size(S,1)) || isinf(size(S,2))
       error("Cannot convert $S to a Matrix")
   end

   if isbanded(S)
       Matrix(BandedMatrix(S))
   elseif isbandedblockbanded(S)
       Matrix(BandedBlockBandedMatrix(S))
   elseif isblockbanded(S)
       Matrix(BlockBandedMatrix(S))
   else
       eltype(S)[S[k,j] for k=1:size(S,1),j=1:size(S,2)]
   end
end

Base.convert(::Type{BandedMatrix},S::Operator) = default_bandedmatrix(S)

function Base.convert(::Type{BlockBandedMatrix},S::Operator)
    if isbandedblockbanded(S)
        BlockBandedMatrix(BandedBlockBandedMatrix(S))
    else
        default_blockbandedmatrix(S)
    end
end

function Base.convert(::Type{RaggedMatrix},S::Operator)
    if isbanded(S)
        RaggedMatrix(BandedMatrix(S))
    elseif isbandedblockbanded(S)
        RaggedMatrix(BandedBlockBandedMatrix(S))
    elseif isblockbanded(S)
        RaggedMatrix(BlockBandedMatrix(S))
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
    elseif isblockbanded(parent(S))
        BlockBandedMatrix(S)
    elseif israggedbelow(parent(S))
        RaggedMatrix(S)
    else
        Matrix(S)
    end
end

Base.convert(::Type{AbstractVector},S::Operator) = Vector(S)
