#### Caching



## Ac'* for QROperatorQ

function mul_pars(Ac::Adjoint{<:Any,<:QROperatorQ},B::AbstractVector,tolerance,maxlength)
    A = parent(Ac)
    T = promote_type(eltype(A),eltype(B))
    mulpars(Operator{T}(A)',Array{T}(B),tolerance,maxlength)
end


include("matrix.jl")
include("ragged.jl")
include("banded.jl")
include("almostbanded.jl")
include("blockbanded.jl")
include("bandedblockbanded.jl")
