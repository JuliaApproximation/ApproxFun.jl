
include("SumSpace.jl")
include("ArraySpace.jl")
include("PiecewiseSpace.jl")
include("ProductSpaceOperators.jl")
include("PrependOperators.jl")
include("SliceSpace.jl")


⊕(A::FunctionSpace,B::FunctionSpace)=domainscompatible(A,B)?SumSpace(A,B):PiecewiseSpace(A,B)
⊕(f::Fun,g::Fun)=Fun(interlace(coefficients(f),coefficients(g)),space(f)⊕space(g))



for TYP in (:SumSpace,:PiecewiseSpace,:FunctionSpace) # Resolve conflict
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
        function coefficients(cfsin::Vector,A::FunctionSpace,B::$TYP)
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
    end
end
