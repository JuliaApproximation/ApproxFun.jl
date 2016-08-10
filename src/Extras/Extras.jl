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

import FastTransforms: pochhammer

include("show.jl")
include("poetry.jl")





include("ReImSpace.jl")
include("simplify.jl")

include("eigs.jl")
include("hacks.jl")
include("fastops.jl")
include("fastcache.jl")

include("autodifferentiation.jl")
include("fractional.jl")


include("dualnumbers.jl")
include("lanczos.jl")
