#### Caching


# back substitution
trtrs!(::Type{Val{'U'}},co::CachedOperator,u::Array) =
                trtrs!(Val{'U'},resizedata!(co,size(u,1),size(u,1)).data,u)




## Ac_mul_B! for QROperatorQ

function Ac_mul_Bpars(A::QROperatorQ,B::Vector,tolerance,maxlength)
    T = promote_type(eltype(A),eltype(B))
    Ac_mul_Bpars(Operator{T}(A),Vector{T}(B),tolerance,maxlength)
end


include("matrix.jl")
include("ragged.jl")
include("banded.jl")
include("almostbanded.jl")
include("bandedblock.jl")
include("bandedblockbanded.jl")
