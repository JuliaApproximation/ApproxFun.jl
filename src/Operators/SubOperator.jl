checkbounds(A::Operator,kr::Colon) = nothing

checkbounds(A::Operator,kr) =
    (maximum(kr) > length(A) || minimum(kr) < 1) && throw(BoundsError(A,kr))


checkbounds(A::Operator,kr::Union{Colon,AbstractCount},jr::Union{Colon,AbstractCount}) = nothing

checkbounds(A::Operator,kr::Union{Colon,AbstractCount},jr) =
    (maximum(jr) > size(A,2) || minimum(jr) < 1) && throw(BoundsError(A,(kr,jr)))

checkbounds(A::Operator,kr,jr::Union{Colon,AbstractCount}) =
    (maximum(kr) > size(A,1)  || minimum(kr) < 1 ) && throw(BoundsError(A,(kr,jr)))

checkbounds(A::Operator,kr,jr) =
    (!isempty(kr) && (maximum(kr) > size(A,1) || minimum(kr) < 1)) ||
    (!isempty(jr) && (maximum(jr) > size(A,2) || minimum(jr) < 1)) &&
    throw(BoundsError(A,(kr,jr)))


checkbounds(A::Operator,K::Block,J::Block) =
     1 ≤ K.n[1] ≤ length(blocklengths(rangespace(A))) && 1 ≤ J.n[1] ≤ length(blocklengths(domainspace(A)))

checkbounds(A::Operator,K::BlockRange,J::BlockRange) =
    isempty(K) || isempty(J) || checkbounds(A,maximum(K),maximum(J))



## SubOperator

struct SubOperator{T,B,I,DI,BI} <: Operator{T}
    parent::B
    indexes::I
    dims::DI
    bandwidths::BI
end



function SubOperator(A,inds,dims,lu)
    checkbounds(A,inds...)
    SubOperator{eltype(A),typeof(A),typeof(inds),
                typeof(dims),typeof(lu)}(A,inds,dims,lu)
end

# work around strange bug with bool size
SubOperator(A,inds,dims::Tuple{Bool,Bool},lu) = SubOperator(A,inds,Int.(dims),lu)

function SubOperator(A,inds::Tuple{Block,Block},lu)
    checkbounds(A,inds...)
    SubOperator(A,inds,(blocklengths(rangespace(A))[inds[1].n[1]],blocklengths(domainspace(A))[inds[2].n[1]]),lu)
end

SubOperator(A,inds::Tuple{Block,Block}) = SubOperator(A,inds,subblockbandwidths(A))
function SubOperator(A,inds::Tuple{<:BlockRange,<:BlockRange})
    checkbounds(A,inds...)
    dims = (sum(blocklengths(rangespace(A))[Int.(inds[1])]),sum(blocklengths(domainspace(A))[Int.(inds[2])]))
    SubOperator(A,inds,dims,(dims[1]-1,dims[2]-1))
end

# cannot infer ranges
SubOperator(A,inds,dims) = SubOperator(A,inds,dims,(dims[1]-1,dims[2]-1))
SubOperator(A,inds) = SubOperator(A,inds,map(length,inds))


convert(::Type{Operator{T}},SO::SubOperator) where {T} =
    SubOperator(Operator{T}(SO.parent),SO.indexes,SO.dims,SO.bandwidths)::Operator{T}

function view(A::Operator,kr::AbstractCount,jr::AbstractCount)
    @assert isinf(size(A,1)) && isinf(size(A,2))
    st=step(kr)
    if isbanded(A) && st==step(jr)  # Otherwise, its not a banded operator
        kr1=first(kr)
        jr1=first(jr)
        l,u=(bandinds(A,1)+kr1-jr1)÷st,(bandinds(A,2)+kr1-jr1)÷st
    else
        l,u=-∞,∞
    end
    SubOperator(A,(kr,jr),size(A),(-l,u))
end


function view(A::Operator,kr::Range,jr::Range)
    st=step(kr)
    if isbanded(A) && st == step(jr)
        kr1=first(kr)
        jr1=first(jr)
        l,u=(bandinds(A,1)+kr1-jr1)÷st,(bandinds(A,2)+kr1-jr1)÷st
        SubOperator(A,(kr,jr),(length(kr),length(jr)),(-l,u))
    else
        SubOperator(A,(kr,jr))
    end
end


function view(A::Operator,kr::UnitRange,jr::UnitRange)
    if isbanded(A)
        shft=first(kr)-first(jr)
        l,u=bandwidth(A,1)-shft,bandinds(A,2)+shft
        SubOperator(A,(kr,jr),(length(kr),length(jr)),(l,u))
    else
        SubOperator(A,(kr,jr))
    end
end

view(A::Operator,::Colon,::Colon) = view(A,1:size(A,1),1:size(A,2))
view(A::Operator,::Colon,jr) = view(A,1:size(A,1),jr)
view(A::Operator,kr,::Colon) = view(A,kr,1:size(A,2))


view(A::Operator,K::Block,J::Block) = SubOperator(A,(K,J))
view(A::Operator,K::Block,j::Colon) = view(A,blockrows(A,K),j)
view(A::Operator,k::Colon,J::Block) = view(A,k,blockcols(A,J))
view(A::Operator,K::Block,j) = view(A,blockrows(A,K),j)
view(A::Operator,k,J::Block) = view(A,k,blockcols(A,J))
view(A::Operator,KR::BlockRange,JR::BlockRange) = SubOperator(A,(KR,JR))

view(A::Operator,k,j) = SubOperator(A,(k,j))



## Needed for Broadcast
Base.Broadcast.containertype(::SubOperator) = Array


reindex(A::Operator, B::Tuple{Block,Any}, kj::Tuple{Any,Any}) =
    (reindex(rangespace(A),(B[1],), (kj[1],))[1], reindex(domainspace(A),tail(B), tail(kj))[1])
# always reindex left-to-right, so if we have only a single tuple, then
# we must be the domainspace
reindex(A::Operator, B::Tuple{Block}, kj::Tuple{Any}) = reindex(domainspace(A),B,kj)

reindex(A::Operator, B::Tuple{AbstractVector{Block},Any}, kj::Tuple{Any,Any}) =
    (reindex(rangespace(A),(B[1],), (kj[1],))[1], reindex(domainspace(A),tail(B), tail(kj))[1])
# always reindex left-to-right, so if we have only a single tuple, then
# we must be the domainspace
reindex(A::Operator, B::Tuple{AbstractVector{Block}}, kj::Tuple{Any}) =
    reindex(domainspace(A),B,kj)
# Blocks are preserved under ranges
for TYP in (:Block,:(AbstractVector{Block}),:(AbstractCount{Block})),
        VTYP in (:AbstractVector,:AbstractCount)
    @eval begin
        reindex(A::Operator, B::Tuple{$VTYP{Int},Any}, kj::Tuple{$TYP,Any}) =
            (reindex(rangespace(A), (B[1],), (kj[1],))[1], reindex(domainspace(A),tail(B), tail(kj))[1])
        reindex(A::Operator, B::Tuple{$VTYP{Int}}, kj::Tuple{$TYP}) =
            reindex(domainspace(A),B,kj)
    end
end



view(A::SubOperator,kr::UnitRange,jr::UnitRange) = view(A.parent,reindex(A,parentindexes(A),(kr,jr))...)
view(A::SubOperator,K::Block,J::Block) = view(A.parent,reindex(A,parentindexes(A),(K,J))...)
function Base.view(A::SubOperator,::Type{FiniteRange},jr::AbstractVector{Int})
    cs = (isbanded(A) || isblockbandedbelow(A)) ? colstop(A,maximum(jr)) : mapreduce(j->colstop(A,j),max,jr)
    view(A,1:cs,jr)
end

view(A::SubOperator,kr,jr) = view(A.parent,reindex(A,parentindexes(A),(kr,jr))...)
view(A::SubOperator,kr::AbstractCount,jr::AbstractCount) = view(A.parent,reindex(A,parentindexes(A),(kr,jr))...)


bandwidth(S::SubOperator,k::Int) = S.bandwidths[k]
bandinds(S::SubOperator) = (-bandwidth(S,1),bandwidth(S,2))
function colstop(S::SubOperator{T,OP,Tuple{UnitRange{Int},UnitRange{Int}}},j::Integer) where {T,OP}
    cs = colstop(parent(S),parentindexes(S)[2][j])
    kr = parentindexes(S)[1]
    n = size(S,1)
    if cs < first(kr)
        0
    elseif cs ≥ last(kr)
        n
    else
        min(n,findfirst(kr,cs))
    end
end
colstart(S::SubOperator{T,OP,Tuple{UnitRange{Int},UnitRange{Int}}},j::Integer) where {T,OP} =
    max(findfirst(parentindexes(S)[1],colstart(parent(S),parentindexes(S)[2][j])),1)
rowstart(S::SubOperator{T,OP,Tuple{UnitRange{Int},UnitRange{Int}}},j::Integer) where {T,OP} =
    max(1,findfirst(parentindexes(S)[2],rowstart(parent(S),parentindexes(S)[1][j])))
rowstop(S::SubOperator{T,OP,Tuple{UnitRange{Int},UnitRange{Int}}},j::Integer) where {T,OP} =
        findfirst(parentindexes(S)[2],rowstop(parent(S),parentindexes(S)[1][j]))


# blocks don't change
blockcolstop(S::SubOperator{T,OP,Tuple{II,JJ}},J::Integer) where {T,OP,II<:Range{Int},JJ<:Range{Int}} =
    blockcolstop(parent(S),J)

israggedbelow(S::SubOperator) = israggedbelow(parent(S))

# since blocks don't change with indexex, neither do blockbandinds
blockbandinds(S::SubOperator{T,OP,Tuple{II,JJ}}) where {T,OP,II<:Range{Int},JJ<:Range{Int}} =
    blockbandinds(parent(S))
function blockbandinds(S::SubOperator{T,B,Tuple{<:BlockRange,<:BlockRange}}) where {T,B}
    KR,JR = parentindexes(S)
    l,u = blockbandinds(parent(S))
    sh = first(KR).n[1]-first(JR).n[1]
    l+sh,u+sh
end

isblockbanded(S::SubOperator{T,B,Tuple{Block,Block}}) where {T,B} = false
isbanded(S::SubOperator{T,B,Tuple{Block,Block}}) where {T,B} = isbandedblockbanded(parent(S))
bandinds(S::SubOperator{T,B,Tuple{Block,Block}}) where {T,B} = subblockbandinds(parent(S))
blockbandinds(S::SubOperator{T,B,Tuple{Block,Block}}) where {T,B} = 0,0

function bbbzeros(S::SubOperator)
    kr,jr=parentindexes(S)
    KO=parent(S)
    l,u=blockbandinds(KO)
    λ,μ=subblockbandinds(KO)

    rt=rangespace(KO)
    dt=domainspace(KO)
    k1,j1=isempty(kr) || isempty(jr) ? (first(kr),first(jr)) :
                                        reindex(S,parentindexes(S),(1,1))

    # each row/column that we differ from the the block start shifts
    # the sub block inds
    J = block(dt,j1)
    K = block(rt,k1)
    jsh=j1-blockstart(dt,J)
    ksh=k1-blockstart(rt,K)

    ret=bbbzeros(eltype(KO),-l,u,-λ+jsh,μ+ksh,
            blocklengths(rangespace(S)),
            blocklengths(domainspace(S)))
end

function bbbzeros(S::SubOperator{T,B,Tuple{<:BlockRange,<:BlockRange}}) where {T,B}
    KR,JR=parentindexes(S)
    KO=parent(S)
    l,u=blockbandinds(KO)::Tuple{Int,Int}
    λ,μ=subblockbandinds(KO)::Tuple{Int,Int}

    rt=rangespace(KO)
    dt=domainspace(KO)
    J=JR[1]
    K=KR[1]
    bl_sh = J.n[1]-K.n[1]

    KBR = blocklengthrange(rt,KR)
    KJR = blocklengthrange(dt,JR)

    ret=bbbzeros(eltype(KO),-l+bl_sh,u-bl_sh,-λ,μ,
            KBR,KJR)
end


function domainspace(S::SubOperator)
    P =parent(S)
    sp=domainspace(P)
    kr=parentindexes(S)[2]

    SubSpace{typeof(sp),typeof(kr),domaintype(sp),rangetype(sp)}(sp,kr)
end
function rangespace(S::SubOperator)
    P =parent(S)
    sp=rangespace(P)
    kr=parentindexes(S)[1]

    SubSpace{typeof(sp),typeof(kr),domaintype(sp),rangetype(sp)}(sp,kr)
end

size(V::SubOperator) = V.dims
size(V::SubOperator,k::Int) = V.dims[k]

unsafe_getindex(V::SubOperator,k::Integer,j::Integer) = V.parent[reindex(V,parentindexes(V),(k,j))...]
getindex(V::SubOperator,k::Integer,j::Integer) = V.parent[reindex(V,parentindexes(V),(k,j))...]
getindex(V::SubOperator,k::Integer,j::Range) = V.parent[reindex(V,parentindexes(V),(k,j))...]
getindex(V::SubOperator,k::Range,j::Integer) = V.parent[reindex(V,parentindexes(V),(k,j))...]
getindex(V::SubOperator,k::Range,j::Range) = V.parent[reindex(V,parentindexes(V),(k,j))...]
Base.parent(S::SubOperator) = S.parent
Base.parentindexes(S::SubOperator) = S.indexes



for OP in (:isblockbanded,:isblockbandedabove,:isblockbandedbelow,
                :isbandedblockbanded,:isbandedblockbandedabove,
                :isbandedblockbandedbelow)
    @eval $OP(S::SubOperator) = $OP(parent(S))
end

# TODO: These should be removed as the general purpose case will work,
# once the notion of bandedness of finite dimensional operators is made sense of


for TYP in (:RaggedMatrix,:Matrix)
    def_TYP = parse("default_" * string(TYP))
    @eval function convert(::Type{$TYP},S::SubOperator)
        if isinf(size(S,1)) || isinf(size(S,2))
            error("Cannot convert $S to a $TYP")
        end

        if isbanded(parent(S))
            $TYP(BandedMatrix(S))
        elseif isbandedblockbanded(parent(S))
            $TYP(BandedBlockBandedMatrix(S))
        elseif isblockbanded(parent(S))
            $TYP(BlockBandedMatrix(S))
        else
            $def_TYP(S)
        end
    end
end

# fast converts to banded matrices would be based on indices, not blocks
function convert(::Type{BandedMatrix},S::SubOperator{T,B,Tuple{<:BlockRange,<:BlockRange}}) where {T,B}
    A = parent(S)
    ds = domainspace(A)
    rs = rangespace(A)
    KR,JR = parentindexes(S)
    BandedMatrix(view(A,
                      blockstart(rs,KR[1]):blockstop(rs,KR[end]),
                      blockstart(ds,JR[1]):blockstop(ds,JR[end])))
end




function A_mul_B_coefficients(A::SubOperator{T,B,Tuple{UnitRange{Int},UnitRange{Int}}},b) where {T,B}
    if size(A,2) == length(b)
        AbstractMatrix(A)*b
    else
        view(A,:,1:length(b))*b
    end
end
