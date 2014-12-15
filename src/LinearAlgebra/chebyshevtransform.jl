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

function plan_chebyshevtransform{T<:Union(Float64,Complex{Float64})}(x::Vector{T})
#TODO confirm that this can handle T=Complex{Float64} looks like this is a real-to-real transform
    if length(x)==1
        return identity
    else
        return FFTW.plan_r2r(x, FFTW.REDFT00)
    end
end

function plan_chebyshevtransform{T<:Number}(x::Vector{T}) return redft00 end

chebyshevtransform(x)=chebyshevtransform(x,plan_chebyshevtransform(x))

function chebyshevtransform{T<:Union(Float64,Complex{Float64})}(x::Vector{T},plan::Function)
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

function chebyshevtransform{T<:Number}(x::Vector{T},plan::Function)
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

#=
chebyshevtransform in terms of ifft, from chebfun
function chebyshevtransform{T<:Number}(x::Vector{T},plan::Function)
    n = length(x)
    if n == 1
        return x
    else
        tmp = cat(1, x[n:-1:2], values[1:n-1]);

        if ( isreal(x) )
            #real-value case
            coeffs = ifft(tmp);
            coeffs = real(coeffs);
        elseif ( isreal(im*x) )
            #Imaginary-valued case:
            coeffs = ifft(imag(tmp));
            coeffs = im*real(coeffs);
        else
            #General case:
            coeffs = ifft(tmp);
        end

        #Truncate:
        coeffs = coeffs[1:n];
        #Scale the interior coefficients:
        coeffs[2:n-1,] = 2*coeffs[2:n-1];
        return coeffs 
    end
end
=#




ichebyshevtransform(x)=ichebyshevtransform(x,plan_chebyshevtransform(x))
function ichebyshevtransform{T<:Number}(x::Vector{T},plan::Function)
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


for func in (:chebyshevtransform,:ichebyshevtransform)
    @eval $func{T<:Integer}(x::Vector{T})=$func(float64(x))
end


function chebyshevtransform(A::Matrix)
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

function ichebyshevtransform(X::Matrix)
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

function chebyshevrootstransform(v::Vector)
    cfs=negateeven!(FFTW.r2r(v,FFTW.REDFT10)/length(v))
    cfs[1]/=2
    cfs    
end

function ichebyshevrootstransform(cfs::Vector)
    cfs[1]*=2
    negateeven!(cfs)
    FFTW.r2r(cfs,FFTW.REDFT01)/2    
end 


