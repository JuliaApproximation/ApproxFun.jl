# This is a Cooley-Tukey FFT algorithm for any number type.
function fft_pow2{T<:Number}(x::Vector{T})
    n = length(x)
    @assert ispow2(n)
    if n==1
        return x
    elseif n==2
        return Complex{real(T)}[x[1]+x[2],x[1]-x[2]]
    end
    even,odd = fft_pow2(x[1:2:end-1]),fft_pow2(x[2:2:end])
    twiddle = exp(-2im*convert(T,π)/n*[0:n-1])
    half1 = even + odd.*twiddle[1:div(n,2)]
    half2 = even + odd.*twiddle[div(n,2)+1:n]
    return vcat(half1,half2)
end
ifft_pow2{T<:Number}(x::Vector{T}) = conj(fft_pow2(conj(x)))/length(x)

# This is a Cooley-Tukey FFT algorithm inspired by many widely available algorithms including:
# c_radix2.c in the GNU Scientific Library and four1 in the Numerical Recipes in C.
# However, the trigonometric recurrence is improved for greater efficiency.
# The algorithm starts with bit-reversal, then divides and conquers in-place.
function fft_pow2!{T<:BigFloat}(x::Vector{T})
    n=length(x)
    nn,j=div(n,2),1
    for i=1:2:n-1
        if j>i
            x[j],x[i]=x[i],x[j]
            x[j+1],x[i+1]=x[i+1],x[j+1]
        end
        m=nn
        while m >= 2 && j > m
            j -= m
            m=div(m,2)
        end
        j += m
    end
    for logn = 2.^collect(1:round(Int,log2(n))-1)
        θ=2convert(T,π)/logn
        wpr,wpi=cos(θ),sin(θ)
        wr,wi=one(T),zero(T)
        for m=1:2:logn-1
            for i=m:2logn:n
                j=i+logn
                mixr,mixi=wr*x[j]-wi*x[j+1],wr*x[j+1]+wi*x[j]
                x[j],x[j+1]=x[i]-mixr,x[i+1]-mixi
                x[i],x[i+1]=x[i]+mixr,x[i+1]+mixi
            end
            wr,wi=wi*wpi+wr*wpr,wi*wpr-wr*wpi
        end
    end
    return x
end

function fft_pow2{T<:BigFloat}(x::Vector{Complex{T}})
    y = interlace(real(x),imag(x))
    fft_pow2!(y)
    return complex(y[1:2:end],y[2:2:end])
end
fft_pow2{T<:BigFloat}(x::Vector{T}) = fft_pow2(complex(x))

function ifft_pow2{T<:BigFloat}(x::Vector{Complex{T}})
    y = interlace(real(x),-imag(x))
    fft_pow2!(y)
    return complex(y[1:2:end],-y[2:2:end])/length(x)
end
