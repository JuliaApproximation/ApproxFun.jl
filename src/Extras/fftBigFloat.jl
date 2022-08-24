# This is a Cooley-Tukey FFT algorithm for any number type.
function fft_pow2(x::StridedVector)
    n = length(x)
    T = mapreduce(eltype(x) isa Fun ? cfstype : eltype, promote_type,x)
    @assert ispow2(n)
    if n==1
        return convert(Vector, x)
    elseif n==2
        return eltype(x)[x[1]+x[2], x[1]-x[2]]
    end
    even,odd = fft_pow2(@view x[1:2:end-1]), fft_pow2(@view x[2:2:end])
    ret = similar(x, n)
    for (halfind, ind) in enumerate(1:div(n,2))
        ret[ind] = even[halfind] + odd[halfind] * cis(-2convert(T,π)/n * (ind - 1))
    end
    for (halfind, ind) in enumerate(div(n,2)+1:n)
        ret[ind] = even[halfind] + odd[halfind] * cis(-2convert(T,π)/n * (ind - 1))
    end
    return ret
end
ifft_pow2(x::StridedVector) = conj(fft_pow2(conj(x)))/length(x)
