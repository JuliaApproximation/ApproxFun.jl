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
SubOperator(A,inds,dims) = SubOperator(A,inds,(dims[1]-1,dims[2]-1))
SubOperator(A,inds) = SubOperator(A,inds,map(length,inds))

function view(A::Operator,kr::AbstractCount,jr::AbstractCount)
    @assert isbanded(A) && isinf(size(A,1)) && isinf(size(A,2))
    st=step(kr)
    @assert st==step(jr)  # Otherwise, its not a banded operator
    kr1=first(kr)
    jr1=first(jr)
    l,u=(bandinds(A,1)+kr1-jr1)÷st,(bandinds(A,2)+kr1-jr1)÷st
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

view(A::Operator,::Colon,jr::Range) = view(A,1:size(A,1),jr)
view(A::Operator,kr::Range,::Colon) = view(A,kr,1:size(A,2))

view(A::SubOperator,kr,jr) =
    view(A.parent,A.indexes[1][kr],A.indexes[2][jr])

bandwidth(S::SubOperator,k::Int) = S.bandwidths[k]
bandinds(S::SubOperator) = (-bandwidth(S,1),bandwidth(S,2))



domainspace(S::SubOperator) = SubSpace(domainspace(S),parentindexes(S)[2])
rangespace(S::SubOperator) = SubSpace(rangespace(S),parentindexes(S)[1])

size(V::SubOperator) = V.dims
unsafe_getindex(V::SubOperator,k::Integer,j::Integer) = V.parent[V.indexes[1][k],V.indexes[2][j]]
getindex(V::SubOperator,k::Integer,j::Integer) = V.parent[V.indexes[1][k],V.indexes[2][j]]
getindex(V::SubOperator,k::Integer,j::Range) = V.parent[V.indexes[1][k],V.indexes[2][j]]
getindex(V::SubOperator,k::Range,j::Integer) = V.parent[V.indexes[1][k],V.indexes[2][j]]
getindex(V::SubOperator,k::Range,j::Range) = V.parent[V.indexes[1][k],V.indexes[2][j]]
Base.parent(S::SubOperator) = S.parent
Base.parentindexes(S::SubOperator) = S.indexes
