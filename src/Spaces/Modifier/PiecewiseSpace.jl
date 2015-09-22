export PiecewiseSpace,depiece,pieces

###########
# Piecewise Space
############

immutable PiecewiseSpace{S,T,d} <: Space{T,AnyDomain,d}
    spaces::Vector{S}
    PiecewiseSpace(::AnyDomain)=new(S[S(AnyDomain())])
    PiecewiseSpace(sp::Vector{S})=new(sp)
end
PiecewiseSpace(sp::Vector{Any})=PiecewiseSpace([sp...])
PiecewiseSpace(sp::Tuple)=PiecewiseSpace([sp...])
PiecewiseSpace(S::Space,spaces::Vector)=PiecewiseSpace{eltype(spaces),basistype(S),ndims(S)}(spaces)
PiecewiseSpace(spaces)=PiecewiseSpace(first(spaces),spaces)
PiecewiseSpace(a::Space,b::Space)=PiecewiseSpace([a,b])




















## space promotion




## cumsum






## Conversion from Vector to Piecewise
#  Right now if a Vector fun has different spaces in each component we represenent
#  it by a piecewise fun, so this allows conversion.

#
# function coefficients(f::Vector,a::VectorSpace,b::PiecewiseSpace)
#     A=a.space
#     n=length(a)
#     @assert n==length(b.spaces)
#     ret=copy(f)
#     for k=1:n
#         ret[k:n:end]=coefficients(ret[k:n:end],A,b.spaces[k])
#     end
#     ret
# end
