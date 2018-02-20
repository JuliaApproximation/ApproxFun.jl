export plan_chebyshevtransform, plan_ichebyshevtransform, chebyshevtransform, ichebyshevtransform

## Transforms take values at Chebyshev points of the first and second kinds and produce Chebyshev coefficients

#TODO confirm that this can handle T=Complex{Float64} looks like this is a real-to-real transform


struct ChebyshevTransformPlan{T,kind,inplace,P} <: Plan{T}
    plan::P
end

ChebyshevTransformPlan{k,inp}(plan) where {k,inp} =
    ChebyshevTransformPlan{eltype(plan),k,inp,typeof(plan)}(plan)



function plan_chebyshevtransform!(x::AbstractVector{T}; kind::Integer=1) where T<:fftwNumber
    if kind == 1
        plan = plan_r2r!(x, REDFT10)
        ChebyshevTransformPlan{1,true}(plan)
    elseif kind == 2
        if length(x) ≤ 1
            error("Cannot create a length $(length(x)) chebyshev transform")
        end
        plan = plan_r2r!(x, REDFT00)
        ChebyshevTransformPlan{2,true}(plan)
    end
end

function plan_chebyshevtransform(x::AbstractVector{T};kind::Integer=1) where T<:fftwNumber
    plan = plan_chebyshevtransform!(x;kind=kind)
    ChebyshevTransformPlan{kind,false}(plan)
end

function *(P::ChebyshevTransformPlan{T,1,true},x::AbstractVector{T}) where T
    n = length(x)
    if n == 1
        x
    else
        x = P.plan*x
        x[1]/=2
        scale!(inv(T(n)),x)
    end
end

function *(P::ChebyshevTransformPlan{T,2,true},x::AbstractVector{T}) where T
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

chebyshevtransform!(x::AbstractVector{T};kind::Integer=1) where {T<:fftwNumber} =
    plan_chebyshevtransform!(x;kind=kind)*x

chebyshevtransform(x;kind::Integer=1) = chebyshevtransform!(copy(x);kind=kind)

*(P::ChebyshevTransformPlan{T,k,false},x::AbstractVector{T}) where {T,k} = P.plan*copy(x)

## Inverse transforms take Chebyshev coefficients and produce values at Chebyshev points of the first and second kinds


struct IChebyshevTransformPlan{T,kind,inplace,P}
    plan::P
end

function plan_ichebyshevtransform!(x::AbstractVector{T};kind::Integer=1) where T<:fftwNumber
    if kind == 1
        if length(x) == 0
            error("Cannot create a length 0 inverse chebyshev transform")
        end
        plan = plan_r2r!(x, REDFT01)
        IChebyshevTransformPlan{T,1,true,typeof(plan)}(plan)
    elseif kind == 2
        if length(x) ≤ 1
            error("Cannot create a length $(length(x)) inverse chebyshev transform")
        end
        plan = plan_chebyshevtransform!(x;kind=2)
        IChebyshevTransformPlan{T,2,true,typeof(plan)}(plan)
    end
end

function plan_ichebyshevtransform(x::AbstractVector{T};kind::Integer=1) where T<:fftwNumber
    plan = plan_ichebyshevtransform!(similar(Vector{T},indices(x));kind=kind)
    IChebyshevTransformPlan{T,kind,false,typeof(plan)}(plan)
end

function *(P::IChebyshevTransformPlan{T,1,true},x::AbstractVector{T}) where T<:fftwNumber
    x[1] *=2
    x = scale!(T(0.5),P.plan*x)
    x
end

function *(P::IChebyshevTransformPlan{T,2,true},x::AbstractVector{T}) where T<:fftwNumber
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

ichebyshevtransform!(x::AbstractVector{T};kind::Integer=1) where {T<:fftwNumber} =
    plan_ichebyshevtransform!(x;kind=kind)*x

ichebyshevtransform(x;kind::Integer=1) = ichebyshevtransform!(copy(x);kind=kind)

*(P::IChebyshevTransformPlan{T,k,false},x::AbstractVector{T}) where {T,k} = P.plan*copy(x)

## Code generation for integer inputs

for func in (:chebyshevtransform,:ichebyshevtransform)
    @eval $func(x::AbstractVector{T};kind::Integer=1) where {T<:Integer} = $func(convert(Float64,x);kind=kind)
end


# Matrix inputs


function chebyshevtransform!(X::AbstractMatrix{T};kind::Integer=1) where T<:fftwNumber
    if kind == 1
        if size(X) == (1,1)
            X
        else
            X=r2r!(X,REDFT10)
            X[:,1]/=2;X[1,:]/=2;
            scale!(1/(size(X,1)*size(X,2)),X)
        end
    elseif kind == 2
        if size(X) == (1,1)
            X
        else
            X=r2r!(X,REDFT00)
            scale!(1/((size(X,1)-1)*(size(X,2)-1)),X)
            X[:,1]/=2;X[:,end]/=2
            X[1,:]/=2;X[end,:]/=2
            X
        end
    end
end

function ichebyshevtransform!(X::AbstractMatrix{T};kind::Integer=1) where T<:fftwNumber
    if kind == 1
        if size(X) == (1,1)
            X
        else
            X[1,:]*=2;X[:,1]*=2
            X = r2r(X,REDFT01)
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
