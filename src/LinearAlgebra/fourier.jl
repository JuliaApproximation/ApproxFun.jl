# diff from CosSpace -> SinSpace

function cosspacediff(v::AbstractVector{T}) where T<:Number
    if length(v)==1
        w = zeros(T,1)
    else
        w = Array{T}(undef, length(v)-1)
        for k=1:length(v)-1
            @inbounds w[k] = -k*v[k+1]
        end
    end

    w
end

# diff from SinSpace -> CosSpace

function sinspacediff(v::AbstractVector{T}) where T<:Number
    w = Array{T}(undef, length(v)+1)
    w[1] = zero(T)
    for k=1:length(v)
        @inbounds w[k+1] = k*v[k]
    end

    w
end

# diff from Fourier -> Fourier

function fourierdiff(v::AbstractVector{T}) where T<:Number
    n = 2(length(v)÷2)+1
    w = Array{T}(undef, n)
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

function taylor_diff(v::AbstractVector{T}) where T<:Number
    w = Array{T}(undef, length(v))
    for k=1:length(v)
        @inbounds w[k] = (k-1)*v[k]
    end

    w
end

# diff from Hardy{false} -> Hardy{false}

function hardyfalse_diff(v::AbstractVector{T}) where T<:Number
    w = Array{T}(undef, length(v))
    for k=1:length(v)
        @inbounds w[k] = -k*v[k]
    end

    w
end

# diff from Laurent -> Laurent

function laurentdiff(v::AbstractVector{T}) where T<:Number
    n = length(v)
    w = Array{T}(undef, n)
    w[1] = zero(T)
    n=length(v)

    for k=1:(isodd(n) ? n÷2 : n÷2-1)
        @inbounds w[2k] = -k*v[2k]
        @inbounds w[2k+1] = k*v[2k+1]
    end

    if iseven(n)
        @inbounds w[n] = -(n÷2)*v[n]
    end

    w
end
