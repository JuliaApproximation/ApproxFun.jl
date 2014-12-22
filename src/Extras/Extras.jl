include("specialfunctions.jl")
include("roots.jl")
include("sample.jl")
include("timeevolution.jl")
include("fftBigFloat.jl")

export digits

#
# This function provides a convenient way to query or specify the BigFloat precision.
# It should be moved to a better file location when BigFloat support is settled.
#
digits(n::Integer) = set_bigfloat_precision(int(ceil(n*log2(10))))
digits() = int(floor(get_bigfloat_precision()*log10(2)))