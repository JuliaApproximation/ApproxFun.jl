export plan_chebyshevtransform, ichebyshevtransform,chebyshevtransform

## transforms


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

function plan_chebyshevtransform{T<:FFTW.fftwNumber}(x::Vector{T})
#TODO confirm that this can handle T=Complex{Float64} looks like this is a real-to-real transform
    if length(x)==1
        return identity
    else
        return FFTW.plan_r2r(x, FFTW.REDFT00)
    end
end

chebyshevtransform{T<:FFTW.fftwNumber}(x::Vector{T})=chebyshevtransform(x,plan_chebyshevtransform(x))

#take values at Chebyshev points and produces Chebyshev coefficients
function chebyshevtransform{T<:FFTW.fftwNumber}(x::Vector{T},plan::Function)
    if length(x) == 1
        x
    else
        ret = plan(x)::typeof(x)
        ret[1] /= 2;ret[end] /= 2   
        negateeven!(ret)
        ret*=1./(length(ret)-1)
        
        ret
    end
end

#following Chebfun's vals2coeffs.m
function chebyshevtransform{T<:Number}(x::Vector{T})
    N = length(x)
    if N== 1
        x
    else
        #has size 2N-1
        tmp = cat(1, reverse(x), x[2:end-1])
        local y
        if T<:Real
            y = real(ifft_gen(tmp))
        else
            y = ifft_gen(tmp)
        end

        y = y[1:N]
        y[2:N-1] *= 2
        return y
    end
end

ichebyshevtransform{T<:FFTW.fftwNumber}(x::Vector{T})=ichebyshevtransform(x,plan_chebyshevtransform(x))
function ichebyshevtransform{T<:FFTW.fftwNumber}(x::Vector{T},plan::Function)
    if(length(x) == 1)
        x
    else
        ##TODO: make thread safe
        x[1] *= 2;x[end] *= 2
        
        ret = chebyshevtransform(x,plan)::typeof(x)
        
        x[1] /=2;x[end] /=2
        
        ret[1] *= 2;ret[end] *= 2
        
        negateeven!(ret)
        
        ret *= .5*(length(x) - 1)
        
        flipud(ret)
    end
end

function ichebyshevtransform{T<:Number}(x::Vector{T})
    if(length(x) == 1)
        x
    else
        ##TODO: make thread safe
        x[1] *= 2;x[end] *= 2
        
        ret = chebyshevtransform(x)::typeof(x)
        
        x[1] /=2;x[end] /=2
        
        ret[1] *= 2;ret[end] *= 2
        
        negateeven!(ret)
        
        ret *= .5*(length(x) - 1)
        
        flipud(ret)
    end
end

function chebyshevtransform{T<:FFTW.fftwNumber}(A::Matrix{T})
    if(size(A) == (1,1))
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
    if(size(X) == (1,1))
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

plan_chebyshevrootstransform(x)=length(x)==1?identity:FFTW.plan_r2r(x, FFTW.REDFT10)

chebyshevrootstransform(x)=chebyshevrootstransform(x,plan_chebyshevrootstransform(x))

function chebyshevrootstransform{T<:FFTW.fftwNumber}(v::Vector{T},plan::Function)
    cfs=negateeven!(plan(v)::typeof(v))
    cfs[1]/=2
    cfs/=length(v)
    cfs    
end

function ichebyshevrootstransform{T<:FFTW.fftwNumber}(cfs::Vector{T})
    cfs[1]*=2
    negateeven!(cfs)
    FFTW.r2r(cfs,FFTW.REDFT01)/2    
end 






for func in (:chebyshevtransform,:ichebyshevtransform,:chebyshevrootstransform,:ichebyshevrootstransform)
    @eval $func{T<:Integer}(x::Vector{T})=$func(float64(x))
end