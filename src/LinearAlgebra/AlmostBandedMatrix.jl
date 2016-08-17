## AlmostBandedMatrix



immutable AlmostBandedMatrix{T} <: AbstractMatrix{T}
    bands::BandedMatrix{T}
    fill::LowRankMatrix{T}
    function AlmostBandedMatrix(bands::BandedMatrix{T},fill::LowRankMatrix{T})
        if size(bands) ≠ size(fill)
            error("Data and fill must be compatible size")
        end
        new(bands,fill)
    end
end

AlmostBandedMatrix(bands::BandedMatrix,fill::LowRankMatrix) =
    AlmostBandedMatrix{promote_type(eltype(bands),eltype(fill))}(bands,fill)

AlmostBandedMatrix{T}(::Type{T},n::Integer,m::Integer,l::Integer,u::Integer,r::Integer) =
    AlmostBandedMatrix(BandedMatrix(T,n,m,l,u),LowRankMatrix(T,n,m,r))


abzeros{T}(::Type{T},n::Integer,m::Integer,l::Integer,u::Integer,r::Integer) =
    AlmostBandedMatrix(bzeros(T,n,m,l,u),lrzeros(T,n,m,r))



Base.size(A::AlmostBandedMatrix) = size(A.bands)
Base.linearindexing{ABM<:AlmostBandedMatrix}(::Type{ABM}) =
    Base.LinearSlow()


function getindex(B::AlmostBandedMatrix,k::Integer,j::Integer)
    if j > k + bandwidth(B.bands,2)
        B.fill[k,j]
    else
        B.bands[k,j]
    end
end

# can only change the bands, not the fill
function Base.setindex!(B::AlmostBandedMatrix,v,k::Integer,j::Integer)
        B.bands[k,j] = v
end


function pad!(B::AlmostBandedMatrix,n::Integer,m::Integer)
    pad!(B.bands,n,m)
    pad!(B.fill,n,m)
    B
end
