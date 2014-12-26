function dft!(CS::Vector,x::Vector)
    N = length(x)
    warn( "Performing a DFT of size $N")
    @assert length(CS) == N
    for k=0:(N-1)
        for n=0:(N-1)
            arg = -im*2pi* k*n/N
            cs = exp(arg)
            @inbounds CS[k+1] += cs*x[n+1]
        end
    end
end

#=  
function Base.fft{T<:Number}(data::Vector{T})
    N = length(data)
    CS = zeros(T,N)
    dft!(CS,data)
    return CS
end
=#

function type_of_real{T<:Number}(x::Type{T})
    if T<:Real
        return T
    else
        return T.parameters[1]
    end
end

function fft_pow2{T<:Number}(x::Vector{T})

    N = length(x)
    @assert ispow2(N)

    if N==1
        return x
    end

    Even = fft_pow2(x[1:2:end-1])
    Odd = fft_pow2(x[2:2:end])

    Tpi = convert(type_of_real(T),pi)

    Twiddle = exp(-2Tpi*im/N * [0:N-1])
    
    half1 = Even + Odd .* Twiddle[1:(N/2)]
    half2 = Even + Odd .* Twiddle[(N/2+1):N]
    return cat(1,half1,half2)
end

#julia 0.3v makes it impossible to overrride fft only for all non-fftw-supported number types.
function fft_gen{T<:Number}(x::Vector{T})
    #=
    if 0==length(T.parameters)
        Tbase = T
    elseif 1==length(T.parameters)
        Tbase = T.parameters[1]
    else
        error("dunno what to do")
    end
    
    S = promote_type(Complex{Tbase},T) #type of output
    
    y = zeros(S,length(x))
    dft!(y,x)
    =#
    N = length(x)
    if ispow2(N)
        return fft_pow2(x)
    end

    #Bluestein's, following http://www.dsprelated.com/dspbooks/mdft/Bluestein_s_FFT_Algorithm.html
    ks = [0:(N-1)]
    ns = [-N:(N-1)]

    Tpi = convert(type_of_real(T),pi)

    W = exp(im*2*Tpi/N)

    Wks = W.^(-(ks.^2)/2)
    xq = x .* Wks
    wq = W.^((ns.^2)/2)

    #convert arrays-to-be-convolved into a common type
    S = promote_type(eltype(xq),eltype(wq))
    xq = convert(Vector{S},xq)
    wq = convert(Vector{S},wq)

    C = conv(xq,wq)[(N+1):(N+N)]

    return Wks .* C
end

function ifft_gen{T<:Number}(x::Vector{T})
    return (1/length(x)) * conj(fft_gen(conj(x)))
end

function ifft_gen!{T<:Number}(x::Vector{T})
    y = (1/length(x)) * conj(fft_gen(conj(x)))
    x[:] = y
    return x
end

#can extend Base.conv appropriately, unlike FFT
function Base.conv{T<:Number}(u::StridedVector{T}, v::StridedVector{T})
    nu = length(u)
    nv = length(v)
    n = nu + nv - 1
    np2 = nextpow2(n)
    upad = [u, zeros(T, np2 - nu)]
    vpad = [v, zeros(T, np2 - nv)]

    Y = fft_gen(upad) .* fft_gen(vpad)
    y = ifft_gen!(Y)
    if T<:Real      #TODO This would not handle Dual/ComplexDual numbers correctly
        y = real(y)
    end

    return y[1:n]
end



function redft00{T<:Number}(x::Vector{T})
    N = length(x)
    warn( "Performing a REDFT00 of size $N because the numbertype is $T")

    y = zeros(T,N,)

    if N==1
        y[1] = x[1]
        return y
    end

    for k=0:(N-1)
        y[k+1] = x[1] + (-1)^k * x[N]

        #will only execute for N>2
        for j=1:(N-2)
            y[k+1] += 2*x[j+1] * cos(pi*j*k/(N-1))
        end
    end
    return y
end 

#= 
b = [BigFloat(r) for r in rand(4000,)]
 abs(fft_gen(b) - fft(b))
=# 
