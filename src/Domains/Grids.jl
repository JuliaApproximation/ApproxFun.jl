

# Represent evenly spaced points on [0,2π)
struct FourierGrid{T} <: AbstractVector{T}
    n::Int
    π_over_n::T
    FourierGrid{T}(n::Int) where T = new{T}(n, convert(T,π)/n)
end

size(g::FourierGrid) = (g.n,)
@inline function getindex(g::FourierGrid, k::Int)
    @boundscheck Base.checkbounds(g, k)
    (2k-2)*g.π_over_n
end
Base.IndexStyle(::Type{<:FourierGrid}) = IndexLinear()


FourierGrid(10)

using BenchmarkTools
g = FourierGrid{Float64}(10_000)

fp(n) = collect(FourierGrid{Float64}(n))
@btime collect(fp(10_000))
@btime fourierpoints(10_000)

n=10
T=Float64
    convert(T,π)*(0:2:2n-2)/n - [convert(T,π)*(2k-2)/n for k=1:n]
