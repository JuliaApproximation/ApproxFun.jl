export plan_chebyshevtransform, plan_ichebyshevtransform, chebyshevtransform, ichebyshevtransform

## Transforms take values at Chebyshev points of the first and second kinds and produce Chebyshev coefficients

#TODO confirm that this can handle T=Complex{Float64} looks like this is a real-to-real transform

function plan_chebyshevtransform{T<:FFTW.fftwNumber}(x::Vector{T};kind::Integer=1)
    if kind == 1
        FFTW.plan_r2r(x, FFTW.REDFT10)
    elseif kind == 2
        if length(x) ≤ 1
            error("Cannot create a length $(length(x)) chebyshev transform")
        end
        FFTW.plan_r2r(x, FFTW.REDFT00)
    end
end

function chebyshevtransform{T<:FFTW.fftwNumber}(x::Vector{T},plan;kind::Integer=1)
    if kind == 1
        n = length(x)
        if n == 1
            x
        else
            ret=negateeven!(plan*x)
            ret[1]/=2
            scale!(inv(T(n)),ret)
        end
    elseif kind == 2
        n = length(x)
        if n == 1
            x
        else
            ret = plan*x
            ret[1] /= 2;ret[end] /= 2
            negateeven!(ret)
            scale!(inv(T(n-1)),ret)
        end
    end
end
chebyshevtransform{T<:FFTW.fftwNumber}(x::Vector{T};kind::Integer=1)=chebyshevtransform(x,plan_chebyshevtransform(x;kind=kind);kind=kind)

## Inverse transforms take Chebyshev coefficients and produce values at Chebyshev points of the first and second kinds


function plan_ichebyshevtransform{T<:FFTW.fftwNumber}(x::Vector{T};kind::Integer=1)
    if kind == 1
        FFTW.plan_r2r(x, FFTW.REDFT01)
    elseif kind == 2
        if length(x) ≤ 1
            error("Cannot create a length $(length(x)) chebyshev transform")
        end

        FFTW.plan_r2r(x, FFTW.REDFT00)
    end
end

function ichebyshevtransform{T<:FFTW.fftwNumber}(x::Vector{T},plan;kind::Integer=1)
    if kind == 1
        x[1] *=2
        ret = scale!(T(0.5),plan*negateeven!(x))
        negateeven!(x)
        x[1]/=2
        ret
    elseif kind == 2
        n = length(x)
        if n == 1
            x
        else
            ##TODO: make thread safe
            x[1] *= 2;x[end] *= 2
            ret = chebyshevtransform(x,plan;kind=kind)
            x[1] /=2;x[end] /=2
            ret[1] *= 2;ret[end] *= 2
            negateeven!(ret)
            scale!(T(.5(n-1)),ret)
            reverse!(ret)
        end
    end
end
ichebyshevtransform{T<:FFTW.fftwNumber}(x::Vector{T};kind::Integer=1)=ichebyshevtransform(x,plan_ichebyshevtransform(x;kind=kind);kind=kind)

## Code generation for integer inputs

for func in (:chebyshevtransform,:ichebyshevtransform)
    @eval $func{T<:Integer}(x::Vector{T};kind::Integer=1)=$func(float64(x);kind=kind)
end


# Matrix inputs


function chebyshevtransform{T<:FFTW.fftwNumber}(X::Matrix{T};kind::Integer=1)
    if kind == 1
        if size(X) == (1,1)
            X
        else
            R=negateeven!(FFTW.r2r(X,FFTW.REDFT10))
            R[:,1]/=2;R[1,:]/=2;
            R/=size(X,1)*size(X,2)
        end
    elseif kind == 2
        if size(X) == (1,1)
            X
        else
            R=FFTW.r2r(X,FFTW.REDFT00)/((size(X,1)-1)*(size(X,2)-1))
            R[:,1]/=2;R[:,end]/=2
            R[1,:]/=2;R[end,:]/=2
            negateeven!(R)
        end
    end
end

function ichebyshevtransform{T<:FFTW.fftwNumber}(X::Matrix{T};kind::Integer=1)
    if kind == 1
        if size(X) == (1,1)
            X
        else
            X[1,:]*=2;X[:,1]*=2
            R = FFTW.r2r(negateeven!(X),FFTW.REDFT01)/4
            negateeven!(X)
            X[1,:]/=2;X[:,1]/=2
            R
        end
    elseif kind == 2
        if size(X) == (1,1)
            X
        else
            X[1,:]*=2;X[end,:]*=2;X[:,1]*=2;X[:,end]*=2
            R=chebyshevtransform(X;kind=kind)
            X[1,:]/=2;X[end,:]/=2;X[:,1]/=2;X[:,end]/=2
            R[1,:]*=2;R[end,:]*=2;R[:,1]*=2;R[:,end]*=2
            negateeven!(R)
            R*=(size(X,1)-1)*(size(X,2)-1)/4
            flipud(fliplr(R))
        end
    end
end
