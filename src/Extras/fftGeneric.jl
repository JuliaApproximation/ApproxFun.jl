typealias BigFloats Union{BigFloat,Complex{BigFloat}}

if VERSION >= v"0.4-dev"
    # old DFT API: p(x) # deprecated
    wrap_fft_plan(x::Function) = x
    # new DFT API
    immutable FFTPlanWrapper{P}
        p::P
    end
    call(p::FFTPlanWrapper, arg) = p.p * arg
    wrap_fft_plan(x) = FFTPlanWrapper(x)
else
    # 0.3 (old) DFT API
    wrap_fft_plan(x) = x
end

# The following implements Bluestein's algorithm, following http://www.dsprelated.com/dspbooks/mdft/Bluestein_s_FFT_Algorithm.html
# To add more types, add them in the union of the function's signature.
function Base.fft{T<:BigFloats}(x::Vector{T})
    n = length(x)
    if ispow2(n) return fft_pow2(x) end
    ks = linspace(zero(real(T)),n-one(real(T)),n)
    Wks = exp(-im*convert(T,π)*ks.^2/n)
    xq,wq = x.*Wks,conj([exp(-im*convert(T,π)*n);reverse(Wks);Wks[2:end]])
    return Wks.*conv(xq,wq)[n+1:2n]
end

Base.ifft{T<:BigFloats}(x::Vector{T}) = conj(fft(conj(x)))/length(x)
function Base.ifft!{T<:BigFloats}(x::Vector{T})
    y = conj(fft(conj(x)))/length(x)
    x[:] = y
    return x
end

function Base.conv{T<:BigFloats}(u::StridedVector{T}, v::StridedVector{T})
    nu,nv = length(u),length(v)
    n = nu + nv - 1
    np2 = nextpow2(n)
    pad!(u,np2),pad!(v,np2)
    y = ifft_pow2(fft_pow2(u).*fft_pow2(v))
    #TODO This would not handle Dual/ComplexDual numbers correctly
    y = T<:Real ? real(y[1:n]) : y[1:n]
end

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

Base.plan_fft{T<:BigFloats}(x::Vector{T}) = fft
Base.plan_ifft{T<:BigFloats}(x::Vector{T}) = ifft
Base.plan_fft{F<:Fun}(x::Vector{F}) = fft
Base.plan_ifft{F<:Fun}(x::Vector{F}) = ifft

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
            ret = ifft([reverse(x),x[2:end-1]])[1:n]
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
            w = exp(-im*convert(T,π)*[0:2n-1]/2n)/2
            w[1] *= 2;w[n+1] *= 0;w[n+2:end] *= -1
            ret = fft(w.*[x,one(T),x[end:-1:2]])[n:-1:1]
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
        v = complex([x[1:n],x[n-1:-1:2]],[0,-x[2n-2:-1:n+1],0,x[n+1:2n-2]])
        real(fft(v))
    end
    plan
end

# SinSpace plans for BigFloat

function plan_transform{T<:BigFloat}(::SinSpace,x::Vector{T})
    function plan(x)
        imag(fft([zero(T),-x,zero(T),reverse(x)]))[2:length(x)+1]
    end
    plan
end

function plan_itransform{T<:BigFloat}(::SinSpace,x::Vector{T})
    function plan(x)
        imag(fft([zero(T),-x,zero(T),reverse(x)]))[2:length(x)+1]
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


### Dual Number support on the spaces SP

for SP in (:Chebyshev,:Fourier,:Laurent,:Taylor,:(Hardy{false}),:SinSpace,:CosSpace)
    @eval begin
        plan_transform{D<:Dual}(S::$SP,vals::Vector{D})=plan_transform(S,value(vals))
        plan_itransform{D<:Dual}(S::$SP,vals::Vector{D})=plan_itransform(S,value(vals))
        transform{D<:Dual}(S::$SP,vals::Vector{D},plan)=Dual(transform(S,value(vals),plan),transform(S,dual(vals),plan))
        itransform{D<:Dual}(S::$SP,vals::Vector{D},plan)=Dual(itransform(S,value(vals),plan),itransform(S,dual(vals),plan))
    end
end
