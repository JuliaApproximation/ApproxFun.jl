
struct SubSpace{DS,IT,DD,RR} <: Space{DD,RR}
    space::DS
    indexes::IT
    SubSpace{DS,IT,DD,RR}(sp::DS,ind::IT) where {DS,IT,DD,RR} = new(sp,ind)
end

SubSpace(sp::Space,kr) =
    SubSpace{typeof(sp),typeof(kr),domaintype(sp),rangetype(sp)}(sp,kr)

SubSpace(sp::SubSpace,kr) = SubSpace(sp.space,reindex(sp,sp.indexes,to_indexes(kr))[1])

domain(DS::SubSpace) = domain(DS.space)
dimension(sp::SubSpace) = length(sp.indexes)

|(sp::Space,kr::AbstractRange) = SubSpace(sp,kr)


function |(f::Fun,kr::AbstractInfUnitRange)
    @assert dimension(space(f)) == ∞
    Fun(space(f)|kr,f.coefficients[kr[1]:end])
end

block(sp::SubSpace, k::Integer) = block(sp.space,reindex(sp,(sp.indexes,),(k,))[1])

function blocklengths(sp::SubSpace{DS,UnitRange{Int}}) where DS
    N = first(sp.indexes)
    M = last(sp.indexes)
    B1=block(sp.space,N)
    B2=block(sp.space,M)
    # if the blocks are equal, we have only one bvlock
    B1 == B2 && return [zeros(Int,B1.n[1]-1);length(sp.indexes)]

    [zeros(Int,B1.n[1]-1);
         blockstop(sp.space,B1)-N+1;blocklengths(sp.space)[B1.n[1]+1:B2.n[1]-1];
        M-blockstart(sp.space,B2)+1]
end

function blocklengths(sp::SubSpace{DS,<:AbstractInfUnitRange{Int}}) where DS
    N = first(sp.indexes)
    B1=block(sp.space,N)

    Vcat([zeros(Int,B1.n[1]-1); blockstop(sp.space,B1)-N+1],
            blocklengths(sp.space)[B1.n[1]+1:∞])
end

blocklengths(sp::SubSpace{DS,Block{1,T}}) where {DS, T} =
    [blocklengths(sp.space)[Int(sp.indexes)]]

blocklengths(sp::SubSpace) = error("Not implemented for non-unitrange subspaces")


## Block reindexing for SubSpace
reindex(sp::SubSpace, b::Tuple{Block{1}}, ks::Tuple{Any}) = (blockstart(sp.space,b[1]).+ks[1].-1,)
reindex(sp::SubSpace, br::Tuple{BlockRange1}, ks::Tuple{Block{1}}) = (Block(Int.(br[1])[first(ks[1].n)]),)
reindex(sp::SubSpace, br::Tuple{BlockRange1}, ks::Tuple{Any}) = (blockstart(sp.space,first(Int.(br[1]))).+ks[1].-1,)

# blocks stay the same with unit range indices
reindex(sp::SubSpace, br::Tuple{AbstractVector{Int}}, ks::Tuple{Block{1}}) =
    reindex(sp, br, (blockrange(sp,first(ks)),))
reindex(sp::SubSpace, br::Tuple{AbstractVector{Int}}, ks::Tuple{BlockRange1}) =
    reindex(sp, br, (blockrange(sp,first(ks)),))


## Block
blocklengths(sp::SubSpace{DS,Block}) where {DS} = [blocklengths(sp.space)[sp.indexes.n[1]]]
dimension(sp::SubSpace{DS,Block}) where {DS} = blocklengths(sp.space)[sp.indexes.n[1]]


blocklengths(sp::SubSpace{DS,BlockRange1}) where {DS} = blocklengths(sp.space)[Int.(sp.indexes)]
dimension(sp::SubSpace{DS,BlockRange1}) where {DS} = sum(blocklengths(sp.space)[Int.(sp.indexes)])


block(_,B::Block) = B



##


spacescompatible(S1::SubSpace{DS,IT,DD,RR},S2::SubSpace{DS,IT,DD,RR}) where {DS,IT,DD,RR} =
    spacescompatible(S1.space,S2.space) && S1.indexes == S2.indexes

==(S1::SubSpace{DS,IT,DD,RR},S2::SubSpace{DS,IT,DD,RR}) where {DS,IT,DD,RR} =
    S1.space == S2.space && S1.indexes == S2.indexes

canonicalspace(a::SubSpace) = a.space



setdomain(DS::SubSpace,d::Domain) = SubSpace(setdomain(DS.space,d),DS.indexes)

Conversion(a::SubSpace,b::Space) = ConcreteConversion(a,b)

function Conversion(a::S,b::SubSpace{S,<:AbstractInfUnitRange{Int}}) where S<:Space
    @assert first(b.indexes) == 1
    ConversionWrapper(SpaceOperator(Operator(I,a),a,b))
end

bandwidths(C::ConcreteConversion{<:SubSpace{<:Any,<:AbstractInfUnitRange{Int}}}) =
    first(domainspace(C).indexes)-1,1-first(domainspace(C).indexes)

getindex(C::ConcreteConversion{<:SubSpace}, k::Integer,j::Integer) =
    domainspace(C).indexes[j]==k ? one(eltype(C)) : zero(eltype(C))


# avoid ambiguity
for OP in (:first,:last)
    @eval getindex(E::ConcreteEvaluation{SubSpace{DS,IT,DD,RR},typeof($OP)},k::Integer) where {IT,DS,DD,RR}=
        Evaluation(E.space.space,E.x,E.order)[E.space.indexes[k]]
end
getindex(E::ConcreteEvaluation{SubSpace{DS,IT,DD,RR}},k::Integer) where {IT,DS,DD,RR}=
    Evaluation(E.space.space,E.x,E.order)[E.space.indexes[k]]
getindex(E::ConcreteEvaluation{SubSpace{DS,IT,DD,RR}},kr::AbstractRange) where {IT,DS,DD,RR} =
    Evaluation(E.space.space,E.x,E.order)[E.space.indexes[kr]]


function conversion_rule(a::SubSpace,b::SubSpace)
     if a==b
        a
    else
        NoSpace()
    end
end
# return the space that has banded Conversion to the other
function conversion_rule(a::SubSpace,b::Space)
    if a.space==b
        a  # we can write droping coefficients as a banded operator
    else
        NoSpace()
    end
end

function union_rule(a::SubSpace,b::SubSpace)
     if a == b
        a
    else
        NoSpace()
    end
end
# return the space that has banded Conversion to the other
function union_rule(a::SubSpace,b::Space)
    if a.space==b
        b  # we can write droping coefficients as a banded operator
    else
        NoSpace()
    end
end


function subspace_coefficients(v::AbstractVector,sp::SubSpace,dropsp::SubSpace)
    if sp == dropsp
        v
    else
        coefficients(v,sp,canonicalspace(sp),dropsp)
    end
end


function subspace_coefficients(v::AbstractVector,sp::Space,dropsp::SubSpace)
    n=length(v)
    if sp == dropsp.space
        ret = Array{eltype(v)}(undef, 0)
        for k in dropsp.indexes
            if k > n
                return ret
            end
            push!(ret,v[k])
        end
        ret
    else
        coefficients(v,sp,canonicalspace(dropsp),dropsp)
    end
end

function subspace_coefficients(v::AbstractVector,dropsp::SubSpace,sp::Space)
    if sp==dropsp.space
        isempty(v) && return v
        ret = zeros(eltype(v),dropsp.indexes[length(v)])
        for k = eachindex(v)
            ret[dropsp.indexes[k]] = v[k]
        end
        ret
    else
        coefficients(v,dropsp,canonicalspace(dropsp),sp)
    end
end


coefficients(v::AbstractVector,sp::SubSpace,dropsp::SubSpace) = subspace_coefficients(v,sp,dropsp)



## points


## transform
function transform(sp::SubSpace, vals::AbstractVector)
    ret = transform(sp.space,vals)
    coefficients(ret,sp.space,sp)
end

itransform(sp::SubSpace,cfs::AbstractVector) =
    itransform(sp.space,coefficients(cfs,sp,sp.space))

for Tran in (:plan_transform, :plan_transform!, :plan_itransform, :plan_itransform!)
    @eval $Tran(sp::SubSpace) = $Tran(sp, dimension(sp))
end

points(sp::SubSpace, n) = points(sp.space, n)
points(sp::SubSpace) = points(sp, dimension(sp))


coefficients(v::AbstractVector,::SubSpace{DS,IT,Segment{Vec{2,TT}}},::TensorSpace{SV,DD}) where {DS,IT,TT,SV,DD<:Domain2d} =
    error("Not callable, only defined for ambiguity errors.")
coefficients(v::AbstractVector,::SubSpace{DS,IT,D},::TensorSpace{SV,DD}) where {DS,IT,D,SV,DD<:Domain2d} =
    error("Not callable, only defined for ambiguity errors.")

for TYP in (:SumSpace,:PiecewiseSpace,:TensorSpace,:ConstantSpace,:Space) # Resolve conflict
    @eval begin
        coefficients(v::AbstractVector,sp::$TYP,dropsp::SubSpace) = subspace_coefficients(v,sp,dropsp)
        coefficients(v::AbstractVector,dropsp::SubSpace,sp::$TYP) = subspace_coefficients(v,dropsp,sp)
    end
end

## ProductFUn

# values{S<:SubSpace,V<:SubSpace}(f::ProductFun{S,V})=values(ProductFun(f,space(f,1).space,space(f,2).space))
# values{S<:SubSpace}(f::ProductFun{S})=values(ProductFun(f,space(f,1).space,space(f,2)))
#
#
# function coefficients{n,DS,TT,DD}(f::ProductFun{SubSpace{n,1,DS,DD,RR}},ox::Space,oy::Space)
#     T=eltype(f)
#     m=size(f,1)
#     A=[pad!(coefficients(fx,ox),m+n) for fx in f.coefficients]
#     B=hcat(A...)::Array{T,2}
#     for k=1:size(B,1)
#         ccfs=coefficients(B[k,:],space(f,2),oy)
#         if length(ccfs)>size(B,2)
#             B=pad(B,size(B,1),length(ccfs))
#         end
#         B[k,1:length(ccfs)]=ccfs
#         #B[k,length(ccfs):1:end]=zero(T)
#     end
#
#     B
# end
