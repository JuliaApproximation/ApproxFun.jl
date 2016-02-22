include("specialfunctions.jl")
include("roots.jl")
include("sample.jl")
include("timeevolution.jl")
include("fftBigFloat.jl")
include("fftGeneric.jl")


#
# This function provides a convenient way to query or specify the BigFloat precision.
#
digits(n::Int) = set_bigfloat_precision(round(Int,ceil(n*log2(10))))
digits() = round(Int,floor(get_bigfloat_precision()*log10(2)))

function pochhammer(x::Number,n::Integer)
    ret = one(x)
    if n≥0
        for i=0:n-1
            ret *= x+i
        end
    else
        ret /= pochhammer(x+n,-n)
    end
    ret
end

pochhammer{T<:Number}(x::AbstractArray{T,1},n::Integer) = [pochhammer(x[i],n) for i=1:length(x)]
pochhammer{T<:Number}(x::AbstractArray{T,2},n::Integer) = [pochhammer(x[i,j],n) for i=1:size(x,1),j=1:size(x,2)]
pochhammer{T<:Number}(x::AbstractArray{T},n::Integer) = reshape([ pochhammer(x[i],n) for i in eachindex(x) ], size(x))

pochhammer(x::Number,n::Number) = gamma(x+n)/gamma(x)
pochhammer{T<:Number}(x::AbstractArray{T},n::Number) = gamma(x+n)./gamma(x)

function pochhammer{T<:Real}(x::Number,n::UnitRange{T})
    ret = Vector{promote_type(typeof(x),T)}(length(n))
    ret[1] = pochhammer(x,first(n))
    for i=2:length(n)
        ret[i] = (x+n[i]-1)*ret[i-1]
    end
    ret
end

include("show.jl")
include("poetry.jl")





include("ReImSpace.jl")
include("simplify.jl")

include("eigs.jl")
include("hacks.jl")

include("autodifferentiation.jl")
include("fractional.jl")
