## SubOperator

immutable SubOperator{T,B,I,DI,BI} <: Operator{T}
    parent::B
    indexes::I
    dims::DI
    bandwidths::BI
end

SubOperator(A,inds,dims,lu) =
    SubOperator{eltype(A),typeof(A),typeof(inds),
                typeof(dims),typeof(lu)}(A,inds,map(length,inds),lu)

# cannot infer ranges
SubOperator(A,inds,dims) = SubOperator(A,inds,dims,(dims[1]-1,dims[2]-1))
SubOperator(A,inds) = SubOperator(A,inds,map(length,inds))

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

view(A::Operator,k,j) = SubOperator(A,(k,j))



view(A::SubOperator,kr::UnitRange,jr::UnitRange) =
    view(A.parent,A.indexes[1][kr],A.indexes[2][jr])



bandwidth(S::SubOperator,k::Int) = S.bandwidths[k]
bandinds(S::SubOperator) = (-bandwidth(S,1),bandwidth(S,2))
function colstop(S::SubOperator,j::Integer)
    cs = colstop(parent(S),parentindexes(S)[2][j])
    kr = parentindexes(S)[2]
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
blockbandinds(S::SubOperator) = (-∞,∞)

function bbbzeros(S::SubOperator)
    kr,jr=parentindexes(S)
    KO=parent(S)
    l,u=blockbandinds(KO)
    λ,μ=subblockbandinds(KO)

    rt=rangetensorizer(KO)
    dt=domaintensorizer(KO)

    J=block(dt,jr[1])
    K=block(rt,kr[1])
    bl_sh = J-K

    jsh=jr[1]-blockstart(dt,J)
    ksh=kr[1]-blockstart(rt,K)
    sbl_sh = jsh-ksh


    ret=bbbzeros(eltype(KO),-l+bl_sh,u-bl_sh,max(-λ,-λ+sbl_sh),max(μ,μ-sbl_sh),
            blocklengthrange(rt,kr),
            blocklengthrange(dt,jr))
end


domainspace(S::SubOperator) = domainspace(parent(S))|parentindexes(S)[2]
rangespace(S::SubOperator) = rangespace(parent(S))|parentindexes(S)[1]

size(V::SubOperator) = V.dims
size(V::SubOperator,k::Int) = V.dims[k]
unsafe_getindex(V::SubOperator,k::Integer,j::Integer) = V.parent[V.indexes[1][k],V.indexes[2][j]]
getindex(V::SubOperator,k::Integer,j::Integer) = V.parent[V.indexes[1][k],V.indexes[2][j]]
getindex(V::SubOperator,k::Integer,j::Range) = V.parent[V.indexes[1][k],V.indexes[2][j]]
getindex(V::SubOperator,k::Range,j::Integer) = V.parent[V.indexes[1][k],V.indexes[2][j]]
getindex(V::SubOperator,k::Range,j::Range) = V.parent[V.indexes[1][k],V.indexes[2][j]]
Base.parent(S::SubOperator) = S.parent
Base.parentindexes(S::SubOperator) = S.indexes




function Base.convert(::Type{RaggedMatrix},S::SubOperator)
    if isbanded(parent(S))
        RaggedMatrix(BandedMatrix(S))
    elseif isbandedblockbanded(parent(S))
        RaggedMatrix(BandedBlockBandedMatrix(S))
    else
        default_raggedmatrix(S)
    end
end
