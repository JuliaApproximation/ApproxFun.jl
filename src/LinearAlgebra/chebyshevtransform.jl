export plan_chebyshevtransform, plan_ichebyshevtransform, chebyshevtransform, ichebyshevtransform

## Transforms take values at Chebyshev points of the first and second kinds and produce Chebyshev coefficients

#TODO confirm that this can handle T=Complex{Float64} looks like this is a real-to-real transform


immutable ChebyshevTransformPlan{kind,T,P}
    plan::P
end



function plan_chebyshevtransform{T<:FFTW.fftwNumber}(x::Vector{T};kind::Integer=1)
    if kind == 1
        plan = FFTW.plan_r2r(x, FFTW.REDFT10)
        ChebyshevTransformPlan{1,T,typeof(plan)}(plan)
    elseif kind == 2
        if length(x) ≤ 1
            error("Cannot create a length $(length(x)) chebyshev transform")
        end
        plan = FFTW.plan_r2r(x, FFTW.REDFT00)
        ChebyshevTransformPlan{2,T,typeof(plan)}(plan)
    end
end

function *{T}(P::ChebyshevTransformPlan{1},x::Vector{T})
    n = length(x)
    if n == 1
        x
    else
        ret = P.plan*x
        ret[1]/=2
        scale!(inv(T(n)),ret)
    end
end

function *{T}(P::ChebyshevTransformPlan{2},x::Vector{T})
    n = length(x)
    if n == 1
        x
    else
        n = length(x)
        if n == 1
            x
        else
            ret = P.plan*x
            ret[1] /= 2;ret[end] /= 2
            scale!(inv(T(n-1)),ret)
        end
    end
end

chebyshevtransform{T<:FFTW.fftwNumber}(x::Vector{T};kind::Integer=1) =
    plan_chebyshevtransform(x;kind=kind)*x

## Inverse transforms take Chebyshev coefficients and produce values at Chebyshev points of the first and second kinds


immutable IChebyshevTransformPlan{kind,T,P}
    plan::P
end

function plan_ichebyshevtransform{T<:FFTW.fftwNumber}(x::Vector{T};kind::Integer=1)
    if kind == 1
        if length(x) == 0
            error("Cannot create a length 0 inverse chebyshev transform")
        end
        plan = FFTW.plan_r2r(x, FFTW.REDFT01)
        IChebyshevTransformPlan{1,T,typeof(plan)}(plan)
    elseif kind == 2
        if length(x) ≤ 1
            error("Cannot create a length $(length(x)) inverse chebyshev transform")
        end
        plan = plan_chebyshevtransform(x;kind=2)
        IChebyshevTransformPlan{2,T,typeof(plan)}(plan)
    end
end

function *{T<:FFTW.fftwNumber}(P::IChebyshevTransformPlan{1},x::Vector{T})
    x[1] *=2
    ret = scale!(T(0.5),P.plan*x)
    x[1]/=2
    ret
end

function *{T<:FFTW.fftwNumber}(P::IChebyshevTransformPlan{2},x::Vector{T})
    n = length(x)
    if n == 1
        x
    else
        ##TODO: make thread safe
        x[1] *= 2;x[end] *= 2
        ret = P.plan*x
        x[1] /=2;x[end] /=2
        ret[1] *= 2;ret[end] *= 2
        scale!(T(.5(n-1)),ret)
    end
end

ichebyshevtransform{T<:FFTW.fftwNumber}(x::Vector{T};kind::Integer=1) =
    plan_ichebyshevtransform(x;kind=kind)*x

## Code generation for integer inputs

for func in (:chebyshevtransform,:ichebyshevtransform)
    @eval $func{T<:Integer}(x::Vector{T};kind::Integer=1) = $func(float64(x);kind=kind)
end


# Matrix inputs


function chebyshevtransform{T<:FFTW.fftwNumber}(X::Matrix{T};kind::Integer=1)
    if kind == 1
        if size(X) == (1,1)
            X
        else
            R=FFTW.r2r(X,FFTW.REDFT10)
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
            R
        end
    end
end

function ichebyshevtransform{T<:FFTW.fftwNumber}(X::Matrix{T};kind::Integer=1)
    if kind == 1
        if size(X) == (1,1)
            X
        else
            X[1,:]*=2;X[:,1]*=2
            R = FFTW.r2r(X,FFTW.REDFT01)/4
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
            R*=(size(X,1)-1)*(size(X,2)-1)/4
            R
        end
    end
end
