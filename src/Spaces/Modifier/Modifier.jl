
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
            #TODO: What if we can convert?  FOr example, A could be Ultraspherical{1}
            # and B could contain Chebyshev
            for k=1:m
                if isconvertible(A,B.spaces[k])
                    cfs=coefficients(cfsin,A,B.spaces[k])
                    ret=zeros(eltype(cfs),m*(length(cfs)-1)+k)
                    ret[k:m:end]=cfs
                    return ret
                end
            end

            defaultcoefficients(cfsin,A,B)
        end
        coefficients(cfsin::Vector,A::Space,B::$TYP)=sumspacecoefficients(cfsin,A,B)
    end
end
