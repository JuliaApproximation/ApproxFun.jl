#### Caching



## Ac_mul_B! for QROperatorQ

function Ac_mul_Bpars(A::QROperatorQ,B::AbstractVector,tolerance,maxlength)
    T = promote_type(eltype(A),eltype(B))
    Ac_mul_Bpars(Operator{T}(A),Array{T}(B),tolerance,maxlength)
end


include("matrix.jl")
include("ragged.jl")
include("banded.jl")
include("almostbanded.jl")
include("blockbanded.jl")
include("bandedblockbanded.jl")
