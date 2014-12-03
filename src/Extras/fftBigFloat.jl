function Base.fft!(data::Vector{BigFloat})
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
        wr=BigFloat("1.0")
        wi=BigFloat("0.0")
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

function Base.fft(x::Vector{Complex{BigFloat}})
    n=length(x)
    y = Array(BigFloat,2n)
    for i = 1:n
        y[2i-1] = x[i].re
        y[2i] = x[i].im
    end
    Base.fft!(y)
    return complex(y[1:2:2n-1],y[2:2:2n])
end
Base.fft(x::Vector{BigFloat}) = Base.fft(complex(x)) #TODO: simplify

function Base.ifft(x::Vector{Complex{BigFloat}})
    n=length(x)
    y = Array(BigFloat,2n)
    for i = 1:n
        y[2i-1] = x[i].re
        y[2i] = -x[i].im
    end
    Base.fft!(y)
    x = complex(y[1:2:2n-1],-y[2:2:2n])/n
end