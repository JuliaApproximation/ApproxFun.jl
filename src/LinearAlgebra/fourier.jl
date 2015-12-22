# diff from CosSpace -> SinSpace

function cosspacediff{T<:Number}(v::Vector{T})
    if length(v)==1
        w = zeros(T,1)
    else
        w = Array{T}(length(v)-1)
        for k=1:length(v)-1
            @inbounds w[k] = -k*v[k+1]
        end
    end

    w
end

# diff from SinSpace -> CosSpace

function sinspacediff{T<:Number}(v::Vector{T})
    w = Array{T}(length(v)+1)
    w[1] = zero(T)
    for k=1:length(v)
        @inbounds w[k+1] = k*v[k]
    end

    w
end

# diff from Fourier -> Fourier

function fourierdiff{T<:Number}(v::Vector{T})
    n = 2(length(v)÷2)+1
    w = Array{T}(n)
    w[1] = zero(T)
    n > 1 && (w[n-1] = zero(T))
    for k=1:n÷2-1
        @inbounds w[2k] = -k*v[2k+1]
        @inbounds w[2k+1] = k*v[2k]
    end
    n > 1 && (w[n] = (n÷2)*v[n-1]; n == length(v) && (w[n-1] = -(n÷2)*v[n]))

    w
end

# diff from Taylor -> Taylor

function taylor_diff{T<:Number}(v::Vector{T})
    w = Array{T}(length(v))
    for k=1:length(v)
        @inbounds w[k] = (k-1)*v[k]
    end

    w
end

# diff from Hardy{false} -> Hardy{false}

function hardyfalse_diff{T<:Number}(v::Vector{T})
    w = Array{T}(length(v))
    for k=1:length(v)
        @inbounds w[k] = -k*v[k]
    end

    w
end

# diff from Laurent -> Laurent

function laurentdiff{T<:Number}(v::Vector{T})
    w = Array{T}(length(v))
    w[1] = zero(T)
    for k=1:length(v)÷2
        @inbounds w[2k] = -k*v[2k]
        @inbounds w[2k+1] = k*v[2k+1]
    end

    w
end
