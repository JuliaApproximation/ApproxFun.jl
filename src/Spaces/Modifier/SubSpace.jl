
immutable SubSpace{DS,IT,T,DD,dim}<: Space{T,DD,dim}
    space::DS
    indexes::IT
end

SubSpace{T,DD,dim}(sp::Space{T,DD,dim},kr) =
    SubSpace{typeof(sp),typeof(kr),T,DD,dim}(sp,kr)

SubSpace(sp::Space,kr) =
    SubSpace{typeof(sp),typeof(kr),basistype(sp),domaintype(sp),dimension(sp)}(sp,kr)

SubSpace(sp::SubSpace,kr) = SubSpace(sp.space,reindex(sp,kr))

domain(DS::SubSpace) = domain(DS.space)
dimension(sp::SubSpace) = length(sp.indexes)
reindex(sp::SubSpace, ks) = sp.indexes[ks]


|(sp::Space,kr::UnitCount) = first(kr)==1?sp:SubSpace(sp,kr)
|(sp::Space,kr::Union{AbstractCount,Range}) = SubSpace(sp,kr)


function |(f::Fun,kr::UnitCount)
    @assert dimension(space(f)) == ∞
    Fun(space(f)|kr,f.coefficients[kr[1]:end])
end

block(sp::SubSpace,k::Integer) = block(sp.space,reindex(sp,k))

function blocklengths{DS}(sp::SubSpace{DS,UnitRange{Int}})
    B1=block(sp.space,sp.indexes[1])
    B2=block(sp.space,sp.indexes[end])
    # if the blocks are equal, we have only one bvlock
    B1 == B2 && return [zeros(Int,B1.K-1);length(sp.indexes)]

    [zeros(Int,B1.K-1);
         blockstop(sp.space,B1)-sp.indexes[1]+1;blocklengths(sp.space)[B1.K+1:B2.K-1];
        sp.indexes[end]-blockstart(sp.space,B2)+1]
end

blocklengths(sp::SubSpace) = error("Not implemented for non-unitrange subspaces")


## Block
blocklengths{DS}(sp::SubSpace{DS,Block}) = [blocklengths(sp.space)[sp.indexes.K]]
reindex{DS}(sp::SubSpace{DS,Block}, ks) = blockstart(sp.space,sp.indexes)+ks-1
dimension{DS}(sp::SubSpace{DS,Block}) = blocklengths(sp.space)[sp.indexes.K]


blocklengths{DS}(sp::SubSpace{DS,UnitRange{Block}}) = blocklengths(sp.space)[Int.(sp.indexes)]
reindex{DS}(sp::SubSpace{DS,UnitRange{Block}}, ks)::Int = blockstart(sp.space,sp.indexes[1])+ks-1

dimension{DS}(sp::SubSpace{DS,UnitRange{Block}}) = sum(blocklengths(sp.space)[Int.(sp.indexes)])


block(::,B::Block) = B
block(::,B::SubBlock) = B.block



##


spacescompatible{DS,IT,T,DD,d}(S1::SubSpace{DS,IT,T,DD,d},S2::SubSpace{DS,IT,T,DD,d}) =
    spacescompatible(S1.space,S2.space) && S1.indexes == S2.indexes

=={DS,IT,T,DD,d}(S1::SubSpace{DS,IT,T,DD,d},S2::SubSpace{DS,IT,T,DD,d}) =
    S1.space == S2.space && S1.indexes == S2.indexes

canonicalspace(a::SubSpace) = a.space



setdomain(DS::SubSpace,d::Domain) = SubSpace(setdomain(DS.space,d),DS.indexes)

Conversion{S<:Space,IT,T,DD,d}(a::SubSpace{S,IT,T,DD,d},b::S) =
    ConcreteConversion(a,b)
bandinds{S,T,DD,d}(C::ConcreteConversion{SubSpace{S,UnitCount{Int},T,DD,d},S}) =
    1-first(domainspace(C).indexes),0

getindex{S,IT,T,DD,d}(C::ConcreteConversion{SubSpace{S,IT,T,DD,d},S},
                   k::Integer,j::Integer) =
    domainspace(C).indexes[j]==k?one(eltype(C)):zero(eltype(C))




getindex{IT,DS,T,DD,d}(E::ConcreteEvaluation{SubSpace{DS,IT,T,DD,d},Bool},k::Integer) =
    Evaluation(E.space.space,E.x,E.order)[E.space.indexes[k]]
getindex{IT,DS,T,DD,d}(E::ConcreteEvaluation{SubSpace{DS,IT,T,DD,d}},k::Integer) =
    Evaluation(E.space.space,E.x,E.order)[E.space.indexes[k]]
getindex{IT,DS,T,DD,d}(E::ConcreteEvaluation{SubSpace{DS,IT,T,DD,d}},kr::Range) =
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


function subspace_coefficients(v::Vector,sp::SubSpace,dropsp::SubSpace)
    if sp == dropsp
        v
    else
        coefficients(v,sp,canonicalspace(sp),dropsp)
    end
end


function subspace_coefficients(v::Vector,sp::Space,dropsp::SubSpace)
    n=length(v)
    if sp == dropsp.space
        ret = Array(eltype(v),0)
        for k in dropsp.indexes
            if k > n
                return ret
            end
            push!(ret,v[k])
        end
    else
        coefficients(v,sp,canonicalspace(dropsp),dropsp)
    end
end

function subspace_coefficients(v::Vector,dropsp::SubSpace,sp::Space)
    if sp==dropsp.space
        ret = zeros(eltype(v),dropsp.indexes[length(v)])
        for k = eachindex(v)
            ret[dropsp.indexes[k]] = v[k]
        end
        ret
    else
        coefficients(v,dropsp,canonicalspace(dropsp),sp)
    end
end


coefficients(v::Vector,sp::SubSpace,dropsp::SubSpace) = subspace_coefficients(v,sp,dropsp)



## points


## transform
function transform(sp::SubSpace,vals::Vector)
    ret=transform(sp.space,vals)
    coefficients(ret,sp.space,sp)
end

itransform(sp::SubSpace,cfs::Vector) =
    itransform(sp.space,coefficients(cfs,sp,sp.space))

points(sp::SubSpace,n) = points(sp.space,n)


coefficients{DS,IT,T,TT,SV,TTT,DD}(v::Vector,::SubSpace{DS,IT,T,Segment{Vec{2,TT}},1},::TensorSpace{SV,TTT,DD,2}) =
    error("Not callable, only defined for ambiguity errors.")
coefficients{DS,IT,T,D,SV,TTT,DD}(v::Vector,::SubSpace{DS,IT,T,D,1},::TensorSpace{SV,TTT,DD,2}) =
    error("Not callable, only defined for ambiguity errors.")

for TYP in (:SumSpace,:PiecewiseSpace,:TensorSpace,:ConstantSpace,:Space) # Resolve conflict
    @eval begin
        coefficients(v::Vector,sp::$TYP,dropsp::SubSpace) = subspace_coefficients(v,sp,dropsp)
        coefficients(v::Vector,dropsp::SubSpace,sp::$TYP) = subspace_coefficients(v,dropsp,sp)
    end
end

## ProductFUn

# values{S<:SubSpace,V<:SubSpace}(f::ProductFun{S,V})=values(ProductFun(f,space(f,1).space,space(f,2).space))
# values{S<:SubSpace}(f::ProductFun{S})=values(ProductFun(f,space(f,1).space,space(f,2)))
#
#
# function coefficients{n,DS,TT,DD}(f::ProductFun{SubSpace{n,1,DS,TT,DD,1}},ox::Space,oy::Space)
#     T=eltype(f)
#     m=size(f,1)
#     A=[pad!(coefficients(fx,ox),m+n) for fx in f.coefficients]
#     B=hcat(A...)::Array{T,2}
#     for k=1:size(B,1)
#         ccfs=coefficients(vec(B[k,:]),space(f,2),oy)
#         if length(ccfs)>size(B,2)
#             B=pad(B,size(B,1),length(ccfs))
#         end
#         B[k,1:length(ccfs)]=ccfs
#         #B[k,length(ccfs):1:end]=zero(T)
#     end
#
#     B
# end
