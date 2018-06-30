export Operator
export bandinds, bandrange, \, periodic
export neumann
export ldirichlet,rdirichlet,lneumann,rneumann
export ldiffbc,rdiffbc,diffbcs
export domainspace,rangespace


abstract type Operator{T} end #T is the entry type, Float64 or Complex{Float64}

eltype(::Operator{T}) where {T} = T
eltype(::Type{Operator{T}}) where {T} = T
eltype(::Type{OT}) where {OT<:Operator} = eltype(supertype(OT))


# default entry type
# we assume entries depend on both the domain and the basis
# realdomain case doesn't use


prectype(sp::Space) = promote_type(prectype(domaintype(sp)),eltype(rangetype(sp)))

 #Operators are struct
copy(A::Operator) = A


## We assume operators are T->T
rangespace(A::Operator) = error("Override rangespace for $(typeof(A))")
domainspace(A::Operator) = error("Override domainspace for $(typeof(A))")
domain(A::Operator) = domain(domainspace(A))


isconstspace(_) = false
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
        Base.size(A::$FF,k::Integer) = k==1 ? 1 : ∞
        ApproxFun.rangespace(F::$FF) = ConstantSpace(eltype(F))
        ApproxFun.isafunctional(::$FF) = true
        ApproxFun.blockbandinds(A::$FF) = 0,hastrivialblocks(domainspace(A)) ? bandinds(A,2) : ∞
        function ApproxFun.defaultgetindex(f::$FF,k::Integer,j::Integer)
            @assert k==1
            f[j]::eltype(f)
        end
        function ApproxFun.defaultgetindex(f::$FF,k::Integer,j::AbstractRange)
            @assert k==1
            f[j]
        end
        function ApproxFun.defaultgetindex(f::$FF,k::Integer,j)
            @assert k==1
            f[j]
        end
        function ApproxFun.defaultgetindex(f::$FF,k::AbstractRange,j::Integer)
            @assert k==1:1
            f[j]
        end
        function ApproxFun.defaultgetindex(f::$FF,k::AbstractRange,j::AbstractRange)
            @assert k==1:1
            reshape(f[j],1,length(j))
        end
        function ApproxFun.defaultgetindex(f::$FF,k::AbstractRange,j)
            @assert k==1:1
            reshape(f[j],1,length(j))
        end
    end
end


nblocks(A::Operator,k) = k==1 ? length(blocklengths(rangespace(A))) : length(blocklengths(domainspace(A)))
nblocks(A::Operator) = (nblocks(A,1),nblocks(A,2))


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


isbandedblockbandedbelow(_) = false
isbandedblockbandedabove(_) = false

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

    return (1-length(blocklengths(rangespace(A))),length(blocklengths(domainspace(A)))-1)
end

# assume dense blocks
subblockbandinds(K::Operator,k) = k==1 ? 1-maximum(blocklengths(rangespace(K))) : maximum(blocklengths(domainspace(K)))-1

isblockbandedbelow(A) = isfinite(blockbandinds(A,1))
isblockbandedabove(A) = isfinite(blockbandinds(A,2))
isblockbanded(A::Operator) = isblockbandedbelow(A) && isblockbandedabove(A)

israggedbelow(A::Operator) = isbandedbelow(A) || isbandedblockbanded(A) || isblockbandedbelow(A)


blockbandwidths(S::Operator) = (-blockbandinds(S,1),blockbandinds(S,2))
blockbandinds(K::Operator,k::Integer) = blockbandinds(K)[k]
blockbandwidth(K::Operator,k::Integer) = k==1 ? -blockbandinds(K,k) : blockbandinds(K,k)

subblockbandwidths(K::Operator) = -subblockbandinds(K,1),subblockbandinds(K,2)
subblockbandinds(K::Operator) = subblockbandinds(K,1),subblockbandinds(K,2)
subblockbandwidth(K::Operator,k::Integer) = k==1 ? -subblockbandinds(K,k) : subblockbandinds(K,k)



bandwidth(A::Operator) = bandwidth(A,1) + bandwidth(A,2) + 1
bandwidth(A::Operator,k::Integer) = k==1 ? -bandinds(A,1) : bandinds(A,2)
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
stride(A::Operator) =
    isdiag(A) ? factorial(10) : 1

isdiag(A::Operator) = bandinds(A)==(0,0)


## Construct operators


include("SubOperator.jl")


#
# sparse(B::Operator,n::Integer)=sparse(BandedMatrix(B,n))
# sparse(B::Operator,n::AbstractRange,m::AbstractRange)=sparse(BandedMatrix(B,n,m))
# sparse(B::Operator,n::Colon,m::AbstractRange)=sparse(BandedMatrix(B,n,m))
# sparse(B::Operator,n::AbstractRange,m::Colon)=sparse(BandedMatrix(B,n,m))

## geteindex



getindex(B::Operator,k,j) = defaultgetindex(B,k,j)
getindex(B::Operator,k) = defaultgetindex(B,k)
getindex(B::Operator,k::Block{2}) = B[Block.(k.n)...]




## override getindex.

defaultgetindex(B::Operator,k::Integer) = error("Override [k] for $(typeof(B))")
defaultgetindex(B::Operator,k::Integer,j::Integer) = error("Override [k,j] for $(typeof(B))")


# Ranges


defaultgetindex(op::Operator,kr::AbstractRange) = eltype(op)[op[k] for k in kr]
defaultgetindex(B::Operator,k::Block,j::Block) = AbstractMatrix(view(B,k,j))
defaultgetindex(B::Operator,k::AbstractRange,j::Block) = AbstractMatrix(view(B,k,j))
defaultgetindex(B::Operator,k::Block,j::AbstractRange) = AbstractMatrix(view(B,k,j))
defaultgetindex(B::Operator,k::AbstractRange,j::AbstractRange) = AbstractMatrix(view(B,k,j))

defaultgetindex(op::Operator,k::Integer,jr::AbstractRange) = eltype(op)[op[k,j] for j in jr]
defaultgetindex(op::Operator,kr::AbstractRange,j::Integer) = eltype(op)[op[k,j] for k in kr]

defaultgetindex(B::Operator,k::Block,j::BlockRange) = AbstractMatrix(view(B,k,j))
defaultgetindex(B::Operator,k::BlockRange,j::BlockRange) = AbstractMatrix(view(B,k,j))

defaultgetindex(op::Operator,k::Integer,jr::BlockRange) = eltype(op)[op[k,j] for j in jr]
defaultgetindex(op::Operator,kr::BlockRange,j::Integer) = eltype(op)[op[k,j] for k in kr]


# Colon casdes
defaultgetindex(A::Operator,kj::CartesianIndex{2}) = A[kj[1],kj[2]]
defaultgetindex(A::Operator,kj::CartesianIndex{1}) = A[kj[1]]
defaultgetindex(A::Operator,k,j) = view(A,k,j)



# TODO: finite dimensional blocks
blockcolstart(A::Operator,J::Integer) = Block(max(1,J-blockbandwidth(A,2)))
blockrowstart(A::Operator,K::Integer) = Block(max(1,K-blockbandwidth(A,1)))
blockcolstop(A::Operator,J::Integer) = Block(min(J+blockbandwidth(A,1),nblocks(A,1)))
blockrowstop(A::Operator,K::Integer) = Block(min(K+blockbandwidth(A,2),nblocks(A,2)))

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
    defOP = Meta.parse("default_"*string(OP))
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

function defaultgetindex(A::Operator,::Type{FiniteRange},jr::BlockRange{1})
    cs = (isbanded(A) || isblockbandedbelow(A)) ? blockcolstop(A,maximum(jr)) : mapreduce(j->blockcolstop(A,j),max,jr)
    A[Block(1):cs,jr]
end

function view(A::Operator,::Type{FiniteRange},jr::AbstractVector{Int})
    cs = (isbanded(A) || isblockbandedbelow(A)) ? colstop(A,maximum(jr)) : mapreduce(j->colstop(A,j),max,jr)
    view(A,1:cs,jr)
end

function view(A::Operator,::Type{FiniteRange},jr::BlockRange{1})
    cs = (isbanded(A) || isblockbandedbelow(A)) ? blockcolstop(A,maximum(jr)) : mapreduce(j->blockcolstop(A,j),max,jr)
    view(A,Block(1):cs,jr)
end


defaultgetindex(A::Operator,K::Block,j::Type{FiniteRange}) = A[blockrows(A,K),j]
defaultgetindex(A::Operator,kr,::Type{FiniteRange}) =
    A[kr,1:rowstop(A,maximum(kr))]





## Composition with a Fun, LowRankFun, and ProductFun

getindex(B::Operator,f::Fun) = B*Multiplication(domainspace(B),f)
getindex(B::Operator,f::LowRankFun{S,M,SS,T}) where {S,M,SS,T} = mapreduce(i->f.A[i]*B[f.B[i]],+,1:rank(f))
getindex(B::Operator{BT},f::ProductFun{S,V,SS,T}) where {BT,S,V,SS,T} =
    mapreduce(i->f.coefficients[i]*B[Fun(f.space[2],[zeros(promote_type(BT,T),i-1);
                                            one(promote_type(BT,T))])],
                +,1:length(f.coefficients))



# Convenience for wrapper ops
unwrap_axpy!(α,P,A) = BLAS.axpy!(α,view(parent(P).op,P.indexes[1],P.indexes[2]),A)
iswrapper(_) = false
haswrapperstructure(_) = false

# use this for wrapper operators that have the same structure but
# not necessarily the same entries
#
#  Ex: c*op or real(op)
macro wrapperstructure(Wrap)
    ret = quote
        haswrapperstructure(::$Wrap) = true
    end

    for func in (:(ApproxFun.bandinds),:(LinearAlgebra.stride),
                 :(ApproxFun.isbandedblockbanded),:(ApproxFun.isblockbanded),
                 :(ApproxFun.israggedbelow),:(Base.size),:(ApproxFun.isbanded),
                 :(ApproxFun.bandwidth),:(ApproxFun.bandwidths),
                 :(ApproxFun.blockbandinds),:(ApproxFun.subblockbandinds),
                 :(LinearAlgebra.issymmetric))
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

        Base.getindex(OP::$Wrap,k::Union{Number,AbstractArray,Colon}...) = OP.op[k...]

        BLAS.axpy!(α,P::ApproxFun.SubOperator{T,OP},A::AbstractMatrix) where {T,OP<:$Wrap} =
            ApproxFun.unwrap_axpy!(α,P,A)

        ApproxFun.mul_coefficients(A::$Wrap,b) = mul_coefficients(A.op,b)
        ApproxFun.mul_coefficients(A::ApproxFun.SubOperator{T,OP,Tuple{UnitRange{Int},UnitRange{Int}}},b) where {T,OP<:$Wrap} =
            mul_coefficients(view(parent(A).op,S.indexes[1],S.indexes[2]),b)
        ApproxFun.mul_coefficients(A::ApproxFun.SubOperator{T,OP},b) where {T,OP<:$Wrap} =
            mul_coefficients(view(parent(A).op,S.indexes[1],S.indexes[2]),b)
    end

    for TYP in (:(BandedMatrices.BandedMatrix),:(ApproxFun.RaggedMatrix),
                :Matrix,:Vector,:AbstractVector)
        ret = quote
            $ret

            Base.convert(::Type{$TYP},P::ApproxFun.SubOperator{T,OP}) where {T,OP<:$Wrap} =
                $TYP(view(parent(P).op,P.indexes[1],P.indexes[2]))
            Base.convert(::Type{$TYP},P::ApproxFun.SubOperator{T,OP,NTuple{2,UnitRange{Int}}}) where {T,OP<:$Wrap} =
                $TYP(view(parent(P).op,P.indexes[1],P.indexes[2]))
        end
    end

    ret = quote
        $ret

        # fast converts to banded matrices would be based on indices, not blocks
        function Base.convert(::Type{BandedMatrices.BandedMatrix},
                              S::ApproxFun.SubOperator{T,OP,NTuple{2,ApproxFun.BlockRange1}}) where {T,OP<:$Wrap}
            A = parent(S)
            ds = domainspace(A)
            rs = rangespace(A)
            KR,JR = parentindexes(S)
            BandedMatrix(view(A,
                              blockstart(rs,first(KR)):blockstop(rs,last(KR)),
                              blockstart(ds,first(JR)):blockstop(ds,last(JR))))
        end


        # if the spaces change, then we need to be smarter
        function Base.convert(::Type{ApproxFun.BlockBandedMatrix},
                              S::ApproxFun.SubOperator{T,OP}) where {T,OP<:$Wrap}
            P = parent(S)
            if blocklengths(domainspace(P)) == blocklengths(domainspace(P.op)) &&
                    blocklengths(rangespace(P)) == blocklengths(rangespace(P.op))
                BlockBandedMatrix(view(parent(S).op,S.indexes[1],S.indexes[2]))
            else
                default_BlockBandedMatrix(S)
            end
        end

        function Base.convert(::Type{ApproxFun.PseudoBlockMatrix},
                              S::ApproxFun.SubOperator{T,OP}) where {T,OP<:$Wrap}
            P = parent(S)
            if blocklengths(domainspace(P)) == blocklengths(domainspace(P.op)) &&
                    blocklengths(rangespace(P)) == blocklengths(rangespace(P.op))
                PseudoBlockMatrix(view(parent(S).op,S.indexes[1],S.indexes[2]))
            else
                default_blockmatrix(S)
            end
        end

        function Base.convert(::Type{ApproxFun.BandedBlockBandedMatrix},
                              S::ApproxFun.SubOperator{T,OP}) where {T,OP<:$Wrap}
            P = parent(S)
            if blocklengths(domainspace(P)) == blocklengths(domainspace(P.op)) &&
                    blocklengths(rangespace(P)) == blocklengths(rangespace(P.op))
                BandedBlockBandedMatrix(view(parent(S).op,S.indexes[1],S.indexes[2]))
            else
                default_BandedBlockBandedMatrix(S)
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



include("ldiv.jl")

include("spacepromotion.jl")
include("banded/banded.jl")
include("general/general.jl")

include("functionals/functionals.jl")
include("almostbanded/almostbanded.jl")

include("systems.jl")

include("qrfact.jl")
include("nullspace.jl")




## Conversion



zero(::Type{Operator{T}}) where {T<:Number} = ZeroOperator(T)
zero(::Type{O}) where {O<:Operator} = ZeroOperator(eltype(O))


eye(S::Space) = IdentityOperator(S)
eye(S::Domain) = eye(Space(S))

convert(A::Type{Operator{T}},f::Fun) where {T} =
    norm(f.coefficients)==0 ? zero(A) : convert(A,Multiplication(f))

convert(A::Type{Operator},f::Fun) =
    norm(f.coefficients)==0 ? ZeroOperator() : Multiplication(f)





## Promotion





promote_rule(::Type{N},::Type{Operator}) where {N<:Number} = Operator{N}
promote_rule(::Type{UniformScaling{N}},::Type{Operator}) where {N<:Number} =
    Operator{N}
promote_rule(::Type{Fun{S,N,VN}},::Type{Operator}) where {S,N<:Number,VN} = Operator{N}
promote_rule(::Type{N},::Type{O}) where {N<:Number,O<:Operator} =
    Operator{promote_type(N,eltype(O))}  # float because numbers are promoted to Fun
promote_rule(::Type{UniformScaling{N}},::Type{O}) where {N<:Number,O<:Operator} =
    Operator{promote_type(N,eltype(O))}
promote_rule(::Type{Fun{S,N,VN}},::Type{O}) where {S,N<:Number,O<:Operator,VN} =
    Operator{promote_type(N,eltype(O))}

promote_rule(::Type{BO1},::Type{BO2}) where {BO1<:Operator,BO2<:Operator} =
    Operator{promote_type(eltype(BO1),eltype(BO2))}




## Wrapper

#TODO: Should cases that modify be included?
const WrapperOperator = Union{SpaceOperator,MultiplicationWrapper,DerivativeWrapper,IntegralWrapper,
                                    ConversionWrapper,ConstantTimesOperator,TransposeOperator}





# The following support converting an Operator to a Matrix or BandedMatrix

## BLAS and matrix routines
# We assume that copy may be overriden

BLAS.axpy!(a, X::Operator, Y::AbstractMatrix) = BLAS.axpy!(a,AbstractMatrix(X),Y)
copyto!(dest::AbstractMatrix, src::Operator) = copyto!(dest, AbstractMatrix(src))

# this is for operators that implement copy via axpy!

BandedMatrix(::Type{Zeros}, V::Operator) = BandedMatrix(Zeros{eltype(V)}(size(V)), bandwidths(V))
Matrix(::Type{Zeros}, V::Operator) = Matrix(Zeros{eltype(V)}(size(V)))
BandedBlockBandedMatrix(::Type{Zeros}, V::Operator) =
    BandedBlockBandedMatrix(Zeros{eltype(V)}(size(V)),
                            (blocklengths(rangespace(V)), blocklengths(domainspace(V))),
                            blockbandwidths(V), subblockbandwidths(V))
BlockBandedMatrix(::Type{Zeros}, V::Operator) =
    BlockBandedMatrix(Zeros{eltype(V)}(size(V)),
                      (AbstractVector{Int}(blocklengths(rangespace(V))),
                       AbstractVector{Int}(blocklengths(domainspace(V)))),
                      blockbandwidths(V))
RaggedMatrix(::Type{Zeros}, V::Operator) =
    RaggedMatrix(Zeros{eltype(V)}(size(V)),
                 Int[max(0,colstop(V,j)) for j=1:size(V,2)])


convert_axpy!(::Type{MT}, S::Operator) where {MT <: AbstractMatrix} =
        BLAS.axpy!(one(eltype(S)), S, MT(Zeros, S))



convert(::Type{BandedMatrix}, S::Operator) = default_BandedMatrix(S)

function convert(::Type{BlockBandedMatrix}, S::Operator)
    if isbandedblockbanded(S)
        BlockBandedMatrix(BandedBlockBandedMatrix(S))
    else
        default_BlockBandedMatrix(S)
    end
end

function default_BlockMatrix(S::Operator)
    ret = PseudoBlockArray(zeros(size(S)),
                        AbstractVector{Int}(blocklengths(rangespace(S))),
                        AbstractVector{Int}(blocklengths(domainspace(S))))
    ret .= S
    ret
end

function convert(::Type{PseudoBlockMatrix}, S::Operator)
    if isbandedblockbanded(S)
        PseudoBlockMatrix(BandedBlockBandedMatrix(S))
    elseif isblockbanded(S)
        PseudoBlockMatrix(BlockBandedMatrix(S))
    else
        default_BlockMatrix(S)
    end
end


# TODO: Unify with SubOperator
for TYP in (:RaggedMatrix, :Matrix)
    def_TYP = Meta.parse("default_" * string(TYP))
    @eval function convert(::Type{$TYP}, S::Operator)
        if isinf(size(S,1)) || isinf(size(S,2))
            error("Cannot convert $S to a $TYP")
        end

        if isbanded(S)
            $TYP(BandedMatrix(S))
        else
            $def_TYP(S)
        end
    end
end

function convert(::Type{Vector}, S::Operator)
    if size(S,2) ≠ 1  || isinf(size(S,1))
        error("Cannot convert $S to a AbstractVector")
    end

    eltype(S)[S[k] for k=1:size(S,1)]
end


# TODO: template out fully
arraytype(::Operator) = Matrix
function arraytype(V::SubOperator{T,B,Tuple{KR,JR}}) where {T, B, KR <: Union{BlockRange, Block}, JR <: Union{BlockRange, Block}}
    P = parent(V)
    isbandedblockbanded(V) && return BandedBlockBandedMatrix
    isblockbanded(V) && return BlockBandedMatrix
    return PseudoBlockMatrix
end

function arraytype(V::SubOperator{T,B,Tuple{KR,JR}}) where {T, B, KR <: Block, JR <: Block}
    P = parent(V)
    isbandedblockbanded(V) && return BandedMatrix
    return Matrix
end


function arraytype(V::SubOperator)
    P = parent(V)
    isbanded(P) && return BandedMatrix
    isinf(size(P,1)) && israggedbelow(P) && return RaggedMatrix
    return Matrix
end

convert(::Type{AbstractMatrix}, V::Operator) = convert(arraytype(V), V)
convert(::Type{AbstractVector}, S::Operator) = Vector(S)
