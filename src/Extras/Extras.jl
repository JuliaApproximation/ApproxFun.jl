
include("sample.jl")

include("fftBigFloat.jl")
include("fftGeneric.jl")


#
# This function provides a convenient way to query or specify the BigFloat precision.
#
digits(n::Int) = set_bigfloat_precision(round(Int,ceil(n*log2(10))))
digits() = round(Int,floor(get_bigfloat_precision()*log10(2)))

include("poetry.jl")

include("simplify.jl")

include("autodifferentiation.jl")
include("fractional.jl")

function dualFun end
function dualcfsFun end

if !isdefined(Base, :get_extension)
	parentdir = dirname(dirname(pathof(@__MODULE__)))
	include(joinpath(parentdir, "ext", "ApproxFunDualNumbersExt.jl"))
end
include("lanczos.jl")
