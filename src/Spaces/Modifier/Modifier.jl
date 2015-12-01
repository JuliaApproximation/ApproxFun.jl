
include("SumSpace.jl")
include("ArraySpace.jl")
include("ProductSpaceOperators.jl")
include("BlockOperators.jl")
include("SliceSpace.jl")


⊕(A::Space,B::Space)=domainscompatible(A,B)?SumSpace(A,B):PiecewiseSpace(A,B)
⊕(f::Fun,g::Fun)=Fun(interlace(coefficients(f),coefficients(g)),space(f)⊕space(g))



# Conversion from Vector to Tuple
# If a Vector fun has different spaces in each component we represenent
#  it by a Tuple fun, so this allows conversion.


function coefficients(f::Vector,a::VectorSpace,b::TupleSpace)
    A=a.space
    n=length(a)
    @assert n==length(b.spaces)
    ret=copy(f)
    for k=1:n
        ret[k:n:end]=coefficients(ret[k:n:end],A,b.spaces[k])
    end
    ret
end





for TYP in (:SumSpace,:PiecewiseSpace,:Space) # Resolve conflict
    @eval begin
        function coefficients(v::Vector,sp::$TYP,dropsp::SliceSpace)
            if sp==dropsp.space
                n=index(dropsp)
                st=stride(dropsp)
                v[st+n:st:end]
            else
                coefficients(v,sp,canonicalspace(dropsp),dropsp)
            end
        end

        function coefficients{V}(v::Vector{V},dropsp::SliceSpace,sp::$TYP)
            if sp==dropsp.space
                n=index(dropsp)
                st=stride(dropsp)
                ret=zeros(V,st*length(v)+n)
                ret[st+n:st:end]=v
                ret
            else
                coefficients(v,dropsp,canonicalspace(dropsp),sp)
            end
        end
    end
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


# add constant to each piece
function coefficients(cfsin::Vector,A::ConstantSpace,B::PiecewiseSpace)
    m=length(B.spaces)

    #assume first that  each space needs 1 coefficient to rep constant
    ret=zeros(eltype(cfsin),m)

    for k=1:m
        cfs=coefficients(cfsin,A,B.spaces[k])
        mm=m*(length(cfs)-1)+k  # max ret, should normally be k
        ret=pad!(ret,max(length(ret),mm))
        ret[k:m:mm]=cfs
    end

    ret
end

for TYP in (:SumSpace,:PiecewiseSpace)
    @eval begin
        # we need to be able to call this to avoid confusion
        function sumspacecoefficients(cfsin::Vector,A::Space,B::$TYP)
            m=length(B.spaces)

            for k=1:m
                if isconvertible(A,B[k])
                    cfs=coefficients(cfsin,A,B[k])
                    ret=zeros(eltype(cfs),k-1)  # all spaces have dimension at least 1,
                                                # so there are k-1 zero coefficients
                                                # corresponding to all spaces less than
                                                # k

                    j=k    # represents the row of ret to be added
                    row=1  # represents the row of cfs to be added

                    @assert length(cfs) ≤ dimension(B[k])
                    dk=length(cfs)

                    # we only allow at most dimension
                    while row ≤ dk
                        push!(ret,cfs[row])  # this sets ret[j] to cfs[row]
                        for λ = k+1:m
                            if dimension(B[λ]) ≥ row # loop through the rest of the spaces
                                push!(ret,0)
                                j+=1
                            end
                        end
                        row += 1   # move on to the next row
                        for λ = 1:k-1
                            if dimension(B[λ]) ≥ row # loop through the previous spaces
                                push!(ret,0)
                                j+=1
                            end
                        end
                        j+=1  # we always increment by 1 for the current space
                    end

                    return ret
                end
            end

            defaultcoefficients(cfsin,A,B)
        end
        coefficients(cfsin::Vector,A::Space,B::$TYP)=sumspacecoefficients(cfsin,A,B)
    end
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
            @assert isa(A[k],Functional)
        end
        @assert isa(A[end],BandedOperator)
        L=LowRankOperator(A[1:end-1],$TYP)
        BB=A[end]
        S=SpaceOperator(StrideOperator(BB,length(A)-1,0),domainspace(BB),
                            $TYP(map(rangespace,A)))
        L+S
    end
end

LowRankPertOperator{OT<:Operator}(A::Vector{OT})=LowRankPertOperator(A,TupleSpace)

function LowRankOperator{FT<:Functional}(Bin::Vector{FT},::Type{PiecewiseSpace})
    B=promotedomainspace(Bin)
    rsp=PiecewiseSpace(map(rangespace,B))
    LowRankOperator(
        Fun{typeof(rsp),Float64}[Fun([zeros(k-1);1],rsp) for k=1:length(B)],
        B)
end

function LowRankOperator{FT<:Functional}(Bin::Vector{FT},::Type{TupleSpace})
    B=promotedomainspace(Bin)
    rsp=TupleSpace(tuple(map(rangespace,B)...,ZeroSpace()))  #TODO: Why the hack?
    LowRankOperator(
        Fun{typeof(rsp),Float64}[Fun([zeros(k-1);1],rsp) for k=1:length(B)],
        B)
end



LowRankOperator{FT<:Functional}(Bin::Vector{FT})=LowRankOperator(Bin,TupleSpace)
