const BigFloats = Union{BigFloat,Complex{BigFloat}}

function fft(x::AbstractVector{F}) where F<:Fun
    n,T = length(x),mapreduce(eltype,promote_type,x)
    if ispow2(n) return fft_pow2(x) end
    ks = range(zero(real(T)),stop=n-one(real(T)),length=n)
    Wks = exp(-im*convert(T,π)*ks.^2/n)
    xq,wq = x.*Wks,conj([exp(-im*convert(T,π)*n);reverse(Wks);Wks[2:end]])
    return Wks.*conv(xq,wq)[n+1:2n]
end

ifft(x::AbstractVector{F}) where {F<:Fun} = conj(fft(conj(x)))/length(x)
function ifft!(x::AbstractVector{F}) where F<:Fun
    y = conj(fft(conj(x)))/length(x)
    x[:] = y
    return x
end

function conv(u::StridedVector{F}, v::StridedVector) where F<:Fun
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

plan_fft(x::Vector{F}) where {F<:Fun} = fft
plan_ifft(x::Vector{F}) where {F<:Fun} = ifft
plan_ifft!(x::Vector{F}) where {F<:Fun} = ifft


# Fourier space plans for BigFloat

plan_transform(sp::Fourier{D,R},x::AbstractVector{T}) where {T<:BigFloat,D,R} =
    TransformPlan{T,typeof(sp),false,Nothing}(sp,nothing)
plan_itransform(sp::Fourier{D,R},x::AbstractVector{T}) where {T<:BigFloat,D,R} =
    ITransformPlan{T,typeof(sp),false,Nothing}(sp,nothing)

function *(::TransformPlan{T,Fourier{D,R},false},x::AbstractVector{T}) where {T<:BigFloat,D,R}
    l = length(x); n = div(l+1,2)
    v = fft(x)
    rmul!(v,convert(T,2)/l)
    v[1] /= 2
    mod(l,2) == 1 ? ApproxFunBase.interlace(real(v[1:n]),-imag(v[2:n])) :
      [ApproxFunBase.interlace(real(v[1:n]),-imag(v[2:n]));-real(v[n+1])/2]
end

function *(::ITransformPlan{T,Fourier{D,R},false},x::AbstractVector{T}) where {T<:BigFloat,D,R}
    l = length(x); n = div(l+1,2)
    # v = complex([x[1:n];x[n-1:-1:2]],[0;-x[2n-2:-1:n+1];0;x[n+1:2n-2]])
    v = mod(l,2) == 1 ?
        complex([2x[1];x[3:2:end];x[end:-2:3]],[0;-x[2:2:end];x[end-1:-2:2]]) :
        complex([2x[1];x[3:2:end-1];-2x[end];x[end-1:-2:3]],[0;-x[2:2:end];x[end-2:-2:2]])
    rmul!(v,convert(T,l)/2)
     real(ifft(v))
end

# SinSpace plans for BigFloat

plan_transform(sp::SinSpace{D,R},x::AbstractVector{T}) where {T<:BigFloat,D,R} =
    TransformPlan{T,typeof(sp),false,Nothing}(sp,nothing)
plan_itransform(sp::SinSpace{D,R},x::AbstractVector{T}) where {T<:BigFloat,D,R} =
    ITransformPlan{T,typeof(sp),false,Nothing}(sp,nothing)


function *(::TransformPlan{T,SinSpace{D,R},false},x::AbstractVector{T}) where {T<:BigFloat,D,R}
    v=imag(fft([0;-x;0;reverse(x)]))[2:length(x)+1]
    rmul!(v,convert(T,1)/(length(x)+1))
    v
end

*(::ITransformPlan{T,SinSpace{D,R},false},x::AbstractVector{T}) where {T<:BigFloat,D,R} =
    imag(fft([0;-x;0;reverse(x)]))[2:length(x)+1]/2


# Fourier space & SinSpace plans for Complex{BigFloat}

for SP in (:Fourier,:SinSpace), (pl,TransPlan) in ((:plan_transform,:TransformPlan),
                                                    (:plan_itransform,:ITransformPlan))
    @eval begin
        $pl(sp::$SP{D,R},x::AbstractVector{T}) where {T<:Complex{BigFloat},D,R} =
                $TransPlan(sp,$pl(sp,Array{T}(undef, length(x))),Val{false})
        *(P::$TransPlan{T,$SP{D,R},false},x::Vector{T}) where {T<:Complex{BigFloat},D,R} =
            complex(P.plan*real(x),P.plan*imag(x))
    end
end
