
include("SumSpace.jl")
include("ArraySpace.jl")
include("ProductSpaceOperators.jl")
include("SubSpace.jl")


⊕(A::Space,B::Space) = domainscompatible(A,B) ? SumSpace(A,B) : PiecewiseSpace(A,B)
⊕(f::Fun,g::Fun) = Fun(space(f) ⊕ space(g), interlace(coefficients(f),coefficients(g)))

+(A::Space,B::Space) = A ⊕ B




#split the cfs into component spaces
function coefficients(cfs::AbstractVector,A::SumSpace,B::SumSpace)
    if spacescompatible(A,B)
        cfs
    else
        mapreduce(f->Fun(f,B),+,components(Fun(A,cfs))).coefficients
    end
end
function coefficients(cfs::AbstractVector,A::PiecewiseSpace,B::PiecewiseSpace)
    if spacescompatible(A,B)
        cfs
    else
        mapreduce(f->Fun(f,B),+,components(Fun(A,cfs))).coefficients
    end
end


# spread a single space into a sum space by placing
# its coefficients depending on k
function interlacewithzeros(cfs::AbstractVector,k,it)
    n = length(cfs)

    ret = Array{eltype(cfs)}(0)
    for (K,j) in it
        if j > n
            break
        elseif K == k
            push!(ret,cfs[j])
        else
            push!(ret,0)
        end
    end

    ret
end

interlacewithzeros(cfs::AbstractVector,k,B::DirectSumSpace) = interlacewithzeros(cfs,k,interlacer(B))

function sumspacecoefficients(cfsin::AbstractVector,A::Space,B::SumSpace)
    m=length(components(B))

    for k=1:m
        if isconvertible(A,component(B,k))
            cfs = coefficients(cfsin,A,component(B,k))
            return interlacewithzeros(cfs,k,B)
        end
    end

    defaultcoefficients(cfsin,A,B)
end

function sumspacecoefficients(cfsin::AbstractVector,A::Space,B::PiecewiseSpace)
    m=length(components(B))

    for k=1:m
        if domain(component(B,k)) == domain(A) && isconvertible(A,component(B,k))
            cfs = coefficients(cfsin,A,component(B,k))
            return interlacewithzeros(cfs,k,B)
        end
    end

    defaultcoefficients(cfsin,A,B)
end

for TYP in (:SumSpace,:PiecewiseSpace), ATYP in (:ConstantSpace,:Space)
    @eval coefficients(cfsin::AbstractVector,A::$ATYP,B::$TYP) = sumspacecoefficients(cfsin,A,B)
end


## LowRank Constructors

# convert a vector of functionals and an operator to a LowRnakPertOperator
# the rangespace is a DirectSumSpace specified by ST of the input rangespaces
# the default is a  ArraySpace, but support is there for PiecewiseSpace
# for bcs
for TYP in (:PiecewiseSpace,:VectorSpace)
    @eval function LowRankPertOperator{OT<:Operator}(A::AbstractVector{OT},::Type{$TYP})
        A=promotedomainspace(A)
        for k=1:length(A)-1
            @assert isafunctional(A[k])
        end
        @assert isbanded(A[end])
        L=LowRankOperator(A[1:end-1],$TYP)
        # add zero functionals to shift down
        BB=[fill(ZeroOperator(domainspace(BB),ConstantSpace()),length(A)-1);A[end]]
        S=InterlaceOperator(BB,domainspace(BB),$TYP(map(rangespace,A)))
        L+S
    end
end

LowRankPertOperator{OT<:Operator}(A::AbstractVector{OT})=LowRankPertOperator(A,VectorSpace)

function LowRankOperator{FT<:Operator}(Bin::AbstractVector{FT},::Type{PiecewiseSpace})
    B=promotedomainspace(Bin)
    rsp=PiecewiseSpace(map(rangespace,B))
    LowRankOperator(
        VFun{typeof(rsp),Float64}[Fun(rsp,[zeros(k-1);1]) for k=1:length(B)],
        B)
end

function LowRankOperator{FT<:Operator}(Bin::AbstractVector{FT},::Type{VectorSpace})
    B=promotedomainspace(Bin)
    rsp=Space([map(rangespace,B);ZeroSpace()])  #TODO: Why the hack?
    LowRankOperator(
        VFun{typeof(rsp),Float64}[Fun(rsp,[zeros(k-1);1]) for k=1:length(B)],
        B)
end



LowRankOperator{FT<:Operator}(Bin::AbstractVector{FT}) = LowRankOperator(Bin,VectorSpace)
