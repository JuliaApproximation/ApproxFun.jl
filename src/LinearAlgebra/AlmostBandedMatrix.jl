## AlmostBandedMatrix



immutable AlmostBandedMatrix{T} <: AbstractMatrix{T}
    bands::BandedMatrix{T}
    fill::LowRankMatrix{T}
    function AlmostBandedMatrix(bands::BandedMatrix{T},fill::LowRankMatrix{T})
        if size(data) ≠ size(fill)
            error("Data and fill must be compatible size")
        end
        new(bands,fill)
    end
end

AlmostBandedMatrix(bands::BandedMatrix,fill::LowRankMatrix) =
    AlmostBandedMatrix{promote_type(eltype(bands),eltype(fill))}(bands,fill)


function getindex(B::AlmostBandedMatrix,k::Integer,j::Integer)
    if j > k + bandwidth(B.data,2)
        B.fill[k,j]
    else
        B.bands[k,j]
    end
end
