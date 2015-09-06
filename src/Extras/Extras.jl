include("specialfunctions.jl")
include("roots.jl")
include("sample.jl")
include("timeevolution.jl")
include("fftBigFloat.jl")
include("fftGeneric.jl")


#
# This function provides a convenient way to query or specify the BigFloat precision.
#
digits(n::Integer) = set_bigfloat_precision(int(ceil(n*log2(10))))
digits() = int(floor(get_bigfloat_precision()*log10(2)))

function pochhammer{T<:Number}(x::T,n::Integer)
    ret = one(T)
    if n>=0
        for i=0:n-1
            ret *= x+i
        end
    else
        ret /= pochhammer(x+n,-n)
    end
    ret
end
function pochhammer{T<:Number}(x::Array{T},n::Integer)
    ret = ones(T,size(x))
    if n>=0
        for i=0:n-1
            ret .*= x+i
        end
    else
        ret ./= pochhammer(x+n,-n)
    end
    ret
end
pochhammer{T1<:Number,T2<:Number}(x::T1,n::T2) = gamma(x+n)/gamma(x)
pochhammer{T1<:Number,T2<:Number}(x::Array{T1},n::T2) = gamma(x+n)./gamma(x)

include("show.jl")
include("poetry.jl")





include("ReImSpace.jl")
include("simplify.jl")

include("eigs.jl")
include("hacks.jl")
