
include("SumSpace.jl")
include("ArraySpace.jl")
include("PiecewiseSpace.jl")
include("ProductSpaceOperators.jl")
include("SliceSpace.jl")



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

for TYP in (:SumSpace,:PiecewiseSpace)
    @eval begin
#split the cfs into component spaces
        coefficients(cfs::Vector,A::$TYP,B::$TYP)=mapreduce(f->Fun(f,B),+,vec(Fun(cfs,A))).coefficients

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
