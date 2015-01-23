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
ifft_pow2{T<:Number}(x::Vector{T}) = conj(fft_gen(conj(x)))/length(x)

# The following are specialized routines for BigFloat data from Numerical Recipes in C.
function fft_pow2!(data::Vector{BigFloat})
    @assert ispow2(length(data))
    nn=int(length(data)/2)
    bigπ=big(π)
    n=nn << 1
    j=1
    for i=1:2:n-1
        if j>i
            data[j],data[i]=data[i],data[j]
            data[j+1],data[i+1]=data[i+1],data[j+1]
        end
        m=nn
        while m >= 2 && j > m
            j -= m
            m >>= 1
        end
        j += m
    end
    mmax=2
    while n > mmax
        istep=mmax << 1
        θ=-2bigπ/mmax
        wtemp=sin(θ/2)
        wpr = -2wtemp*wtemp
        wpi=sin(θ)
        wr=big(1.0)
        wi=big(0.0)
        for m=1:2:mmax-1
            for i=m:istep:n
                j=i+mmax
                tempr=wr*data[j]-wi*data[j+1]
                tempi=wr*data[j+1]+wi*data[j]
                data[j]=data[i]-tempr
                data[j+1]=data[i+1]-tempi
                data[i] += tempr
                data[i+1] += tempi
            end
            wr=(wtemp=wr)*wpr-wi*wpi+wr
            wi=wi*wpr+wtemp*wpi+wi
        end
        mmax=istep
    end
    return data
end

function fft_pow2(x::Vector{Complex{BigFloat}})
    @assert ispow2(length(x))
    n=length(x)
    y = interlace(real(x),imag(x))
    fft_pow2!(y)
    return complex(y[1:2:end],y[2:2:end])
end
fft_pow2(x::Vector{BigFloat}) = fft_pow2(complex(x))

function ifft_pow2(x::Vector{Complex{BigFloat}})
    @assert ispow2(length(x))
    n=length(x)
    y = interlace(real(x),-imag(x))
    fft_pow2!(y)
    x = complex(y[1:2:2n-1],-y[2:2:2n])/n
end
