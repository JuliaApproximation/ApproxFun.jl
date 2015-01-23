# The following implements Bluestein's algorithm, following http://www.dsprelated.com/dspbooks/mdft/Bluestein_s_FFT_Algorithm.html
# To add more types, add them in the union of the function's signature.
function Base.fft{T<:Union(BigFloat,Complex{BigFloat})}(x::Vector{T})
    n = length(x)
    if ispow2(n) return fft_pow2(x) end
    ks = linspace(zero(real(T)),n-one(real(T)),n)
    Wks = exp(-im*convert(T,π)*ks.^2/n)
    xq,wq = x.*Wks,conj([exp(-im*convert(T,π)*n),reverse(Wks),Wks[2:end]])
    return Wks.*conv(xq,wq)[n+1:2n]
end

function Base.ifft{T<:Union(BigFloat,Complex{BigFloat})}(x::Vector{T})
    return conj(fft(conj(x)))/length(x)
end

function Base.ifft!{T<:Union(BigFloat,Complex{BigFloat})}(x::Vector{T})
    y = conj(fft(conj(x)))/length(x)
    x[:] = y
    return x
end

function Base.conv{T<:Number}(u::StridedVector{T}, v::StridedVector{T})
    nu,nv = length(u),length(v)
    n = nu + nv - 1
    np2 = nextpow2(n)
    pad!(u,np2),pad!(v,np2)
    y = ifft_pow2(fft_pow2(u).*fft_pow2(v))
    #TODO This would not handle Dual/ComplexDual numbers correctly
    y = T<:Real ? real(y[1:n]) : y[1:n]
end
