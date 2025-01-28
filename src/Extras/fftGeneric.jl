typealias BigFloats Union{BigFloat,Complex{BigFloat}}

function Base.fft{F<:Fun}(x::Vector{F})
    n,T = length(x),mapreduce(eltype,promote_type,x)
    if ispow2(n) return fft_pow2(x) end
    ks = linspace(zero(real(T)),n-one(real(T)),n)
    Wks = exp(-im*convert(T,π)*ks.^2/n)
    xq,wq = x.*Wks,conj([exp(-im*convert(T,π)*n);reverse(Wks);Wks[2:end]])
    return Wks.*conv(xq,wq)[n+1:2n]
end

Base.ifft{F<:Fun}(x::Vector{F}) = conj(fft(conj(x)))/length(x)
function Base.ifft!{F<:Fun}(x::Vector{F})
    y = conj(fft(conj(x)))/length(x)
    x[:] = y
    return x
end

function Base.conv{F<:Fun}(u::StridedVector{F}, v::StridedVector)
    nu,nv = length(u),length(v)
    n = nu + nv - 1
    np2 = nextpow2(n)
    pad!(u,np2),pad!(v,np2)
    y = ifft_pow2(fft_pow2(u).*fft_pow2(v))
    #TODO This would not handle Dual/ComplexDual numbers correctly
    T = promote_type(mapreduce(eltype,promote_type,u),mapreduce(eltype,promote_type,v))
    y = T<:Real ? real(y[1:n]) : y[1:n]
end

######################################################################
# TO BE DEPRECATED FOR V0.4 UPGRADE
######################################################################

# plan_fft for BigFloats (covers Laurent svfft)

Base.plan_fft{F<:Fun}(x::Vector{F}) = fft
Base.plan_ifft{F<:Fun}(x::Vector{F}) = ifft
Base.plan_ifft!{F<:Fun}(x::Vector{F}) = ifft

# Chebyshev transforms and plans for BigFloats

plan_chebyshevtransform{T<:BigFloats}(x::Vector{T};kwds...) = identity
plan_ichebyshevtransform{T<:BigFloats}(x::Vector{T};kwds...) = identity

#following Chebfun's @Chebtech1/vals2coeffs.m and @Chebtech2/vals2coeffs.m
function chebyshevtransform{T<:BigFloats}(x::Vector{T},plan;kind::Integer=1)
    if kind == 1
        n = length(x)
        if n == 1
            x
        else
            w = [2exp(im*convert(T,π)*k/2n) for k=0:n-1]
            ret = w.*ifft([reverse(x);x])[1:n]
            ret = T<:Real ? real(ret) : ret
            ret[1] /= 2
            ret
        end
    elseif kind == 2
        n = length(x)
        if n == 1
            x
        else
            ret = ifft([reverse(x);x[2:end-1]])[1:n]
            ret = T<:Real ? real(ret) : ret
            ret[2:n-1] *= 2
            ret
        end
    end
end

#following Chebfun's @Chebtech1/vals2coeffs.m and @Chebtech2/vals2coeffs.m
function ichebyshevtransform{T<:BigFloats}(x::Vector{T},plan;kind::Integer=1)
    if kind == 1
        n = length(x)
        if n == 1
            x
        else
            w = [exp(-im*convert(T,π)*k/2n)/2 for k=0:2n-1]
            w[1] *= 2;w[n+1] *= 0;w[n+2:end] *= -1
            ret = fft(w.*[x;one(T);x[end:-1:2]])[n:-1:1]
            ret = T<:Real ? real(ret) : ret
        end
    elseif kind == 2
        n = length(x)
        if n == 1
            x
        else
            ##TODO: make thread safe
            x[1] *= 2;x[end] *= 2
            ret = chebyshevtransform(x;kind=kind)
            x[1] /=2;x[end] /=2
            ret[1] *= 2;ret[end] *= 2
            negateeven!(ret)
            ret *= .5*(n-1)
            reverse!(ret)
        end
    end
end

# Fourier space plans for BigFloat

function plan_transform{T<:BigFloat,D}(::Fourier{D},x::Vector{T})
    function plan(x)
        v = fft(x)
        n = div(length(x),2)+1
        [real(v[1:n]);imag(v[n-1:-1:2])]
    end
    plan
end

function plan_itransform{T<:BigFloat,D}(::Fourier{D},x::Vector{T})
    function plan(x)
        n = div(length(x),2)+1
        v = complex([x[1:n];x[n-1:-1:2]],[0;-x[2n-2:-1:n+1];0;x[n+1:2n-2]])
        real(fft(v))
    end
    plan
end

# SinSpace plans for BigFloat

function plan_transform{T<:BigFloat}(::SinSpace,x::Vector{T})
    function plan(x)
        imag(fft([0;-x;0;reverse(x)]))[2:length(x)+1]
    end
    plan
end

function plan_itransform{T<:BigFloat}(::SinSpace,x::Vector{T})
    function plan(x)
        imag(fft([0;-x;0;reverse(x)]))[2:length(x)+1]
    end
    plan
end

# Fourier space & SinSpace plans for Complex{BigFloat}

for SP in (:Fourier,:SinSpace), pl in (:plan_transform,:plan_itransform)
    @eval begin
        function $pl{T<:Complex{BigFloat},D}(::$SP{D},x::Vector{T})
            function plan(x)
                complex($pl($SP(),real(x))(real(x)),$pl($SP(),imag(x))(imag(x)))
            end
            plan
        end
    end
end
