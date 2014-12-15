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

#=      #somehow this method definition is more specific and gets called over built-in fft.
function Base.fft{T<:Number}(data::Vector{T})
    N = length(data)
    CS = zeros(T,N)
    dft!(CS,data)
    return CS
end
=#

function redft00{T<:Number}(x::Vector{T})
    N = length(x)
    warn( "Performing a REDFT00 of size $N")

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



