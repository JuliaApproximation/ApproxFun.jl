
include("SumSpace.jl")
include("ArraySpace.jl")
include("ProductSpaceOperators.jl")
include("BlockOperators.jl")
include("SubSpace.jl")


⊕(A::Space,B::Space)=domainscompatible(A,B)?SumSpace(A,B):PiecewiseSpace(A,B)
⊕(f::Fun,g::Fun)=Fun(interlace(coefficients(f),coefficients(g)),space(f)⊕space(g))



# Conversion from Vector to Tuple
# If a Vector fun has different spaces in each component we represenent
#  it by a Tuple fun, so this allows conversion.


function coefficients(f::Vector,a::VectorSpace,b::TupleSpace)
    error("Reimplement")
    A=a.space
    n=length(a)
    @assert n==length(b.spaces)
    ret=copy(f)
    for k=1:n
        ret[k:n:end]=coefficients(ret[k:n:end],A,b.spaces[k])
    end
    ret
end



#split the cfs into component spaces
function coefficients(cfs::Vector,A::SumSpace,B::SumSpace)
    if spacescompatible(A,B)
        cfs
    else
        mapreduce(f->Fun(f,B),+,vec(Fun(cfs,A))).coefficients
    end
end
function coefficients(cfs::Vector,A::PiecewiseSpace,B::PiecewiseSpace)
    if spacescompatible(A,B)
        cfs
    else
        mapreduce(f->Fun(f,B),+,pieces(Fun(cfs,A))).coefficients
    end
end



function sumspacecoefficients(cfsin::Vector,A::Space,B::SumSpace)
    m=length(B.spaces)

    for k=1:m
        if isconvertible(A,B[k])
            cfs = coefficients(cfsin,A,B[k])

            return interlace([[zeros(B[j]) for j=1:k-1];Fun(cfs,B[k]);[zeros(B[j]) for j=k+1:length(B.spaces)]],
                                B)
        end
    end

    defaultcoefficients(cfsin,A,B)
end

function sumspacecoefficients(cfsin::Vector,A::Space,B::PiecewiseSpace)
    m=length(B.spaces)

    for k=1:m
        if domain(B[k]) == domain(A) && isconvertible(A,B[k])
            cfs = coefficients(cfsin,A,B[k])

            return interlace([[zeros(B[j]) for j=1:k-1];Fun(cfs,B[k]);[zeros(B[j]) for j=k+1:length(B.spaces)]],
                                B)
        end
    end

    defaultcoefficients(cfsin,A,B)
end

for TYP in (:SumSpace,:PiecewiseSpace)
    @eval coefficients(cfsin::Vector,A::Space,B::$TYP) = sumspacecoefficients(cfsin,A,B)
end


## LowRank Constructors

# convert a vector of functionals and an operator to a LowRnakPertOperator
# the rangespace is a DirectSumSpace specified by ST of the input rangespaces
# the default is a  TupleSpace, but support is there for PiecewiseSpace
# for bcs
for TYP in (:PiecewiseSpace,:TupleSpace)
    @eval function LowRankPertOperator{OT<:Operator}(A::Vector{OT},::Type{$TYP})
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

LowRankPertOperator{OT<:Operator}(A::Vector{OT})=LowRankPertOperator(A,TupleSpace)

function LowRankOperator{FT<:Operator}(Bin::Vector{FT},::Type{PiecewiseSpace})
    B=promotedomainspace(Bin)
    rsp=PiecewiseSpace(map(rangespace,B))
    LowRankOperator(
        Fun{typeof(rsp),Float64}[Fun([zeros(k-1);1],rsp) for k=1:length(B)],
        B)
end

function LowRankOperator{FT<:Operator}(Bin::Vector{FT},::Type{TupleSpace})
    B=promotedomainspace(Bin)
    rsp=TupleSpace(tuple(map(rangespace,B)...,ZeroSpace()))  #TODO: Why the hack?
    LowRankOperator(
        Fun{typeof(rsp),Float64}[Fun([zeros(k-1);1],rsp) for k=1:length(B)],
        B)
end



LowRankOperator{FT<:Operator}(Bin::Vector{FT})=LowRankOperator(Bin,TupleSpace)
