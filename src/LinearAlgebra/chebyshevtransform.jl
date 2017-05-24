export plan_chebyshevtransform, plan_ichebyshevtransform, chebyshevtransform, ichebyshevtransform

## Transforms take values at Chebyshev points of the first and second kinds and produce Chebyshev coefficients

#TODO confirm that this can handle T=Complex{Float64} looks like this is a real-to-real transform


struct ChebyshevTransformPlan{T,kind,inplace,P} <: FFTW.Plan{T}
    plan::P
end

(::Type{ChebyshevTransformPlan{k,inp}}){k,inp}(plan) =
    ChebyshevTransformPlan{eltype(plan),k,inp,typeof(plan)}(plan)



function plan_chebyshevtransform!{T<:FFTW.fftwNumber}(x::Vector{T};kind::Integer=1)
    if kind == 1
        plan = FFTW.plan_r2r!(x, FFTW.REDFT10)
        ChebyshevTransformPlan{1,true}(plan)
    elseif kind == 2
        if length(x) ≤ 1
            error("Cannot create a length $(length(x)) chebyshev transform")
        end
        plan = FFTW.plan_r2r!(x, FFTW.REDFT00)
        ChebyshevTransformPlan{2,true}(plan)
    end
end

function plan_chebyshevtransform{T<:FFTW.fftwNumber}(x::Vector{T};kind::Integer=1)
    plan = plan_chebyshevtransform!(x;kind=kind)
    ChebyshevTransformPlan{kind,false}(plan)
end

function *{T}(P::ChebyshevTransformPlan{T,1,true},x::Vector{T})
    n = length(x)
    if n == 1
        x
    else
        x = P.plan*x
        x[1]/=2
        scale!(inv(T(n)),x)
    end
end

function *{T}(P::ChebyshevTransformPlan{T,2,true},x::Vector{T})
    n = length(x)
    if n == 1
        x
    else
        n = length(x)
        if n == 1
            x
        else
            x = P.plan*x
            x[1] /= 2;x[end] /= 2
            scale!(inv(T(n-1)),x)
        end
    end
end

chebyshevtransform!{T<:FFTW.fftwNumber}(x::Vector{T};kind::Integer=1) =
    plan_chebyshevtransform!(x;kind=kind)*x

chebyshevtransform(x;kind::Integer=1) = chebyshevtransform!(copy(x);kind=kind)

*{T,k}(P::ChebyshevTransformPlan{T,k,false},x::Vector{T}) = P.plan*copy(x)

## Inverse transforms take Chebyshev coefficients and produce values at Chebyshev points of the first and second kinds


struct IChebyshevTransformPlan{T,kind,inplace,P}
    plan::P
end

function plan_ichebyshevtransform!{T<:FFTW.fftwNumber}(x::Vector{T};kind::Integer=1)
    if kind == 1
        if length(x) == 0
            error("Cannot create a length 0 inverse chebyshev transform")
        end
        plan = FFTW.plan_r2r!(x, FFTW.REDFT01)
        IChebyshevTransformPlan{T,1,true,typeof(plan)}(plan)
    elseif kind == 2
        if length(x) ≤ 1
            error("Cannot create a length $(length(x)) inverse chebyshev transform")
        end
        plan = plan_chebyshevtransform!(x;kind=2)
        IChebyshevTransformPlan{T,2,true,typeof(plan)}(plan)
    end
end

function plan_ichebyshevtransform{T<:FFTW.fftwNumber}(x::Vector{T};kind::Integer=1)
    plan = plan_ichebyshevtransform!(x;kind=kind)
    IChebyshevTransformPlan{T,kind,false,typeof(plan)}(plan)
end

function *{T<:FFTW.fftwNumber}(P::IChebyshevTransformPlan{T,1,true},x::Vector{T})
    x[1] *=2
    x = scale!(T(0.5),P.plan*x)
    x
end

function *{T<:FFTW.fftwNumber}(P::IChebyshevTransformPlan{T,2,true},x::Vector{T})
    n = length(x)
    if n == 1
        x
    else
        ##TODO: make thread safe
        x[1] *= 2;x[end] *= 2
        x = P.plan*x
        x[1] *= 2;x[end] *= 2
        scale!(T(.5(n-1)),x)
    end
end

ichebyshevtransform!{T<:FFTW.fftwNumber}(x::Vector{T};kind::Integer=1) =
    plan_ichebyshevtransform!(x;kind=kind)*x

ichebyshevtransform(x;kind::Integer=1) = ichebyshevtransform!(copy(x);kind=kind)

*{T,k}(P::IChebyshevTransformPlan{T,k,false},x::Vector{T}) = P.plan*copy(x)

## Code generation for integer inputs

for func in (:chebyshevtransform,:ichebyshevtransform)
    @eval $func{T<:Integer}(x::Vector{T};kind::Integer=1) = $func(convert(Float64,x);kind=kind)
end


# Matrix inputs


function chebyshevtransform!{T<:FFTW.fftwNumber}(X::Matrix{T};kind::Integer=1)
    if kind == 1
        if size(X) == (1,1)
            X
        else
            X=FFTW.r2r!(X,FFTW.REDFT10)
            X[:,1]/=2;X[1,:]/=2;
            scale!(1/(size(X,1)*size(X,2)),X)
        end
    elseif kind == 2
        if size(X) == (1,1)
            X
        else
            X=FFTW.r2r!(X,FFTW.REDFT00)
            scale!(1/((size(X,1)-1)*(size(X,2)-1)),X)
            X[:,1]/=2;X[:,end]/=2
            X[1,:]/=2;X[end,:]/=2
            X
        end
    end
end

function ichebyshevtransform!{T<:FFTW.fftwNumber}(X::Matrix{T};kind::Integer=1)
    if kind == 1
        if size(X) == (1,1)
            X
        else
            X[1,:]*=2;X[:,1]*=2
            X = FFTW.r2r(X,FFTW.REDFT01)
            scale!(1/4,X)
        end
    elseif kind == 2
        if size(X) == (1,1)
            X
        else
            X[1,:]*=2;X[end,:]*=2;X[:,1]*=2;X[:,end]*=2
            X=chebyshevtransform!(X;kind=kind)
            X[1,:]*=2;X[end,:]*=2;X[:,1]*=2;X[:,end]*=2
            scale!((size(X,1)-1)*(size(X,2)-1)/4,X)
        end
    end
end
