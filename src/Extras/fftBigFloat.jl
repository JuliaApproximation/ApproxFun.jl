# This is a Cooley-Tukey FFT algorithm for any number type.
function fft_pow2{F}(x::Vector{F})
    n = length(x)
    T = mapreduce(eltype,promote_type,x)
    @assert ispow2(n)
    if n==1
        return x
    elseif n==2
        return F[x[1]+x[2];x[1]-x[2]]
    end
    even,odd = fft_pow2(x[1:2:end-1]),fft_pow2(x[2:2:end])
    twiddle = exp(-2im*convert(T,π)/n*collect(0:n-1))
    half1 = even + odd.*twiddle[1:div(n,2)]
    half2 = even + odd.*twiddle[div(n,2)+1:n]
    return vcat(half1,half2)
end
ifft_pow2(x::Vector) = conj(fft_pow2(conj(x)))/length(x)
