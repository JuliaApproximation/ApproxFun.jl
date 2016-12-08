checkbounds(A::Operator,kr::Colon) = nothing

checkbounds(A::Operator,kr) =
    (maximum(kr) > length(A) || minimum(kr) < 1) && throw(BoundsError(A,kr))


checkbounds(A::Operator,kr::Union{Colon,AbstractCount},jr::Union{Colon,AbstractCount}) = nothing

checkbounds(A::Operator,kr::Union{Colon,AbstractCount},jr) =
    (maximum(jr) > size(A,2) || minimum(jr) < 1) && throw(BoundsError(A,(kr,jr)))

checkbounds(A::Operator,kr,jr::Union{Colon,AbstractCount}) =
    (maximum(kr) > size(A,1)  || minimum(kr) < 1 ) && throw(BoundsError(A,(kr,jr)))

checkbounds(A::Operator,kr,jr) =
    (maximum(kr) > size(A,1) || maximum(jr) > size(A,2) ||
     minimum(kr) < 1 || minimum(jr) < 1) && throw(BoundsError(A,(kr,jr)))


checkbounds(A::Operator,K::Block,J::Block) =
     1 ≤ K.K ≤ length(blocklengths(rangespace(A))) && 1 ≤ J.K ≤ length(blocklengths(domainspace(A)))

checkbounds(A::Operator,K::Range{Block},J::Range{Block}) =
     checkbounds(A,maximum(K),maximum(J))



## SubOperator

immutable SubOperator{T,B,I,DI,BI} <: Operator{T}
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

function SubOperator(A,inds::Tuple{Block,Block},lu)
    checkbounds(A,inds...)
    SubOperator(A,inds,(blocklengths(rangespace(A))[inds[1].K],blocklengths(domainspace(A))[inds[2].K]),lu)
end

SubOperator(A,inds::Tuple{Block,Block}) = SubOperator(A,inds,subblockbandwidths(A))
function SubOperator(A,inds::Tuple{UnitRange{Block},UnitRange{Block}})
    checkbounds(A,inds...)
    dims = (sum(blocklengths(rangespace(A))[Int.(inds[1])]),sum(blocklengths(domainspace(A))[Int.(inds[2])]))
    SubOperator(A,inds,dims,(dims[1]-1,dims[2]-1))
end

# cannot infer ranges
SubOperator(A,inds,dims) = SubOperator(A,inds,dims,(dims[1]-1,dims[2]-1))
SubOperator(A,inds) = SubOperator(A,inds,map(length,inds))


Base.convert{T}(::Type{Operator{T}},SO::SubOperator) =
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
        l,u=max(bandwidth(A,1)-shft,0),max(bandinds(A,2)+shft,0)
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
view(A::Operator,KR::UnitRange{Block},JR::UnitRange{Block}) = SubOperator(A,(KR,JR))

view(A::Operator,k,j) = SubOperator(A,(k,j))



view(A::SubOperator,kr::UnitRange,jr::UnitRange) =
    view(A.parent,A.indexes[1][kr],A.indexes[2][jr])



bandwidth(S::SubOperator,k::Int) = S.bandwidths[k]
bandinds(S::SubOperator) = (-bandwidth(S,1),bandwidth(S,2))
function colstop(S::SubOperator,j::Integer)
    cs = colstop(parent(S),parentindexes(S)[2][j])
    kr = parentindexes(S)[1]
    n = size(S,1)
    if cs < first(kr)
        1
    elseif cs ≥ last(kr)
        n
    else
        min(n,findfirst(kr,cs))
    end
end
colstart(S::SubOperator,j::Integer) =
    max(findfirst(parentindexes(S)[1],colstart(parent(S),parentindexes(S)[2][j])),1)
rowstart(S::SubOperator,j::Integer) =
    max(1,findfirst(parentindexes(S)[2],rowstart(parent(S),parentindexes(S)[1][j])))
rowstop(S::SubOperator,j::Integer) =
        findfirst(parentindexes(S)[2],rowstop(parent(S),parentindexes(S)[1][j]))


# blocks don't change
blockcolstop(S::SubOperator,J::Integer) = blockcolstop(parent(S),J)

israggedbelow(S::SubOperator) = israggedbelow(parent(S))

# since blocks don't change, neither do blockbandinds
blockbandinds(S::SubOperator) = blockbandinds(parent(S))

function bbbzeros(S::SubOperator)
    kr,jr=parentindexes(S)
    KO=parent(S)
    l,u=blockbandinds(KO)
    λ,μ=subblockbandinds(KO)

    rt=rangetensorizer(KO)
    dt=domaintensorizer(KO)

    k1,j1=reindex(S,(1,1))
    J=block(dt,j1)
    K=block(rt,k1)
    bl_sh = J.K-K.K

    # each row/column that we differ from the the block start shifts
    # the sub block inds
    jsh=j1-blockstart(dt,J)
    ksh=k1-blockstart(rt,K)

    ret=bbbzeros(eltype(KO),-l+bl_sh,u-bl_sh,-λ+jsh,μ+ksh,
            blocklengthrange(rt,kr),
            blocklengthrange(dt,jr))
end


function domainspace(S::SubOperator)
    P =parent(S)
    sp=domainspace(P)
    kr=parentindexes(S)[2]

    SubSpace{typeof(sp),typeof(kr),basistype(sp),typeof(domain(P)),dimension(sp)}(sp,kr)
end
function rangespace(S::SubOperator)
    P =parent(S)
    sp=rangespace(P)
    kr=parentindexes(S)[1]

    SubSpace{typeof(sp),typeof(kr),basistype(sp),typeof(domain(P)),dimension(sp)}(sp,kr)
end

size(V::SubOperator) = V.dims
size(V::SubOperator,k::Int) = V.dims[k]
reindex(V::SubOperator,kj::Tuple{Int,Int}) =
    (reindex(rangespace(V),kj[1])::Int,reindex(domainspace(V),kj[2])::Int)
reindex(V::SubOperator,kj) = (reindex(rangespace(V),kj[1]),reindex(domainspace(V),kj[2]))
    
unsafe_getindex(V::SubOperator,k::Integer,j::Integer) = V.parent[reindex(V,(k,j))...]
getindex(V::SubOperator,k::Integer,j::Integer) = V.parent[reindex(V,(k,j))...]
getindex(V::SubOperator,k::Integer,j::Range) = V.parent[reindex(V,(k,j))...]
getindex(V::SubOperator,k::Range,j::Integer) = V.parent[reindex(V,(k,j))...]
getindex(V::SubOperator,k::Range,j::Range) = V.parent[reindex(V,(k,j))...]
Base.parent(S::SubOperator) = S.parent
Base.parentindexes(S::SubOperator) = S.indexes



for OP in (:isbandedblock,:isbandedblockabove,:isbandedblockbelow,
                :isbandedblockbanded,:isbandedblockbandedabove,
                :isbandedblockbelow)
    @eval $OP(S::SubOperator) = $OP(parent(S))
end

# TODO: These should be removed as the general purpose case will work,
# once the notion of bandedness of finite dimensional operators is made sense of
function Base.convert(::Type{RaggedMatrix},S::SubOperator)
    if isbanded(parent(S))
        RaggedMatrix(BandedMatrix(S))
    elseif isbandedblockbanded(parent(S))
        RaggedMatrix(BandedBlockBandedMatrix(S))
    elseif isbandedblock(parent(S))
        RaggedMatrix(BandedBlockMatrix(S))
    else
        default_raggedmatrix(S)
    end
end
