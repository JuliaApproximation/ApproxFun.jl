export plan_chebyshevtransform, chebyshevtransform, ichebyshevtransform
export plan_chebyshevrootstransform, chebyshevrootstransform, ichebyshevrootstransform

## transforms



#TODO confirm that this can handle T=Complex{Float64} looks like this is a real-to-real transform
plan_chebyshevtransform{T<:FFTW.fftwNumber}(x::Vector{T})=length(x)==1?identity:FFTW.plan_r2r(x, FFTW.REDFT00)

chebyshevtransform{T<:FFTW.fftwNumber}(x::Vector{T})=chebyshevtransform(x,plan_chebyshevtransform(x))

#take values at Chebyshev points and produces Chebyshev coefficients
function chebyshevtransform{T<:FFTW.fftwNumber}(x::Vector{T},plan::Function)
    n = length(x)
    if n == 1
        x
    else
        ret = plan(x)::typeof(x)
        ret[1] /= 2;ret[end] /= 2
        negateeven!(ret)./(n-1)
    end
end

#following Chebfun's @Chebtech2/vals2coeffs.m
function chebyshevtransform{T<:Number}(x::Vector{T})
    n = length(x)
    if n == 1
        x
    else
        ret = ifft([reverse(x),x[2:end-1]])[1:n]
        ret = T<:Real ? real(ret) : ret
        ret[2:n-1] *= 2
        ret
    end
end

ichebyshevtransform{T<:FFTW.fftwNumber}(x::Vector{T})=ichebyshevtransform(x,plan_chebyshevtransform(x))
function ichebyshevtransform{T<:FFTW.fftwNumber}(x::Vector{T},plan::Function)
    n = length(x)
    if n == 1
        x
    else
        ##TODO: make thread safe
        x[1] *= 2;x[end] *= 2
        ret = chebyshevtransform(x,plan)::typeof(x)
        x[1] /=2;x[end] /=2
        ret[1] *= 2;ret[end] *= 2
        negateeven!(ret)
        ret *= .5*(n-1)
        reverse!(ret)
    end
end

function ichebyshevtransform{T<:Number}(x::Vector{T})
    n = length(x)
    if n == 1
        x
    else
        ##TODO: make thread safe
        x[1] *= 2;x[end] *= 2
        ret = chebyshevtransform(x)::typeof(x)
        x[1] /=2;x[end] /=2
        ret[1] *= 2;ret[end] *= 2
        negateeven!(ret)
        ret *= .5*(n-1)
        reverse!(ret)
    end
end

function chebyshevtransform{T<:FFTW.fftwNumber}(A::Matrix{T})
    if size(A) == (1,1)
        A
    else
        R=FFTW.r2r(A,FFTW.REDFT00)/((size(A,1)-1)*(size(A,2)-1))
        R[:,1]/=2;R[:,end]/=2
        R[1,:]/=2;R[end,:]/=2
        negateeven!(R)
        R
    end
end

function ichebyshevtransform{T<:Number}(X::Matrix{T})
    if size(X) == (1,1)
        X
    else
        X[1,:]*=2;X[end,:]*=2;X[:,1]*=2;X[:,end]*=2
        R=chebyshevtransform(X)
        X[1,:]/=2;X[end,:]/=2;X[:,1]/=2;X[:,end]/=2
        R[1,:]*=2;R[end,:]*=2;R[:,1]*=2;R[:,end]*=2
        negateeven!(R)
        R*=(size(X,1)-1)*(size(X,2)-1)/4
        flipud(fliplr(R))
    end
end

## First kind transform


plan_chebyshevrootstransform{T<:FFTW.fftwNumber}(x::Vector{T})=length(x)==1?identity:FFTW.plan_r2r(x, FFTW.REDFT10)

chebyshevrootstransform{T<:FFTW.fftwNumber}(x::Vector{T})=chebyshevrootstransform(x,plan_chebyshevrootstransform(x))

function chebyshevrootstransform{T<:FFTW.fftwNumber}(x::Vector{T},plan::Function)
    ret=negateeven!(plan(x)::typeof(x))
    ret[1]/=2
    ret/=length(x)
end

function ichebyshevrootstransform{T<:FFTW.fftwNumber}(x::Vector{T})
    ret = deepcopy(x)
    ret[1] *=2
    negateeven!(ret)
    FFTW.r2r(ret,FFTW.REDFT01)/2
end

#following Chebfun's @Chebtech1/vals2coeffs.m
function chebyshevrootstransform{T<:Number}(x::Vector{T})
    n = length(x)
    if n == 1
        x
    else
        w = 2exp(im*convert(T,π)*[0:n-1]/2n)
        ret = w.*ifft([reverse(x),x])[1:n]
        ret = T<:Real ? real(ret) : ret
        ret[1] /= 2
        ret
    end
end

function ichebyshevrootstransform{T<:Number}(x::Vector{T})
    n = length(x)
    if n == 1
        x
    else
        w = exp(-im*convert(T,π)*[0:2n-1]/2n)/2
        w[1] *= 2;w[n+1] *= 0;w[n+2:end] *= -1
        ret = fft(w.*[x,one(T),x[end:-1:2]])[n:-1:1]
        ret = T<:Real ? real(ret) : ret
    end
end


for func in (:chebyshevtransform,:ichebyshevtransform,:chebyshevrootstransform,:ichebyshevrootstransform)
    @eval $func{T<:Integer}(x::Vector{T})=$func(float64(x))
end


# Helper routines

function negateeven!(x::Vector)
    for k =2:2:length(x)
        x[k] = -x[k]
    end

    x
end

#checkerboard, same as applying negativeeven! to all rows then all columns
function negateeven!(X::Matrix)
    for k =2:2:size(X,1),j=1:2:size(X,2)
        X[k,j] *= -1
    end
    for k =1:2:size(X,1),j=2:2:size(X,2)
        X[k,j] *= -1
    end

    X
end
