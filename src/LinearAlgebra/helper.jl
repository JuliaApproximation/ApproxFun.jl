export dotu
import Base.chop

# Used for spaces not defined yet
immutable UnsetNumber <: Number  end
Base.promote_rule{N<:Number}(::Type{UnsetNumber},::Type{N})=N

# Test the number of arguments a function takes
if VERSION < v"0.5.0-dev"
    hasnumargs(f,k)=(isgeneric(f)&&applicable(f,zeros(k)...)) || (!isgeneric(f)&&arglength(f)==k)
else
    hasnumargs(f,k)=applicable(f,zeros(k)...)
end


isapprox(a...;kwds...)=Base.isapprox(a...;kwds...)
isapprox(a::Vec,b::Vec;kwds...)=isapprox([a...],[b...];kwds...)

# This creates ApproxFun.real, ApproxFun.eps and ApproxFun.dou
# which we override for default julia types
real(x...)=Base.real(x...)
real(::Type{UnsetNumber})=UnsetNumber
real{T<:Real}(::Type{T})=T
real{T<:Real}(::Type{Complex{T}})=T
real{T<:Real,n}(::Type{Array{T,n}})=Array{T,n}
real{T<:Complex,n}(::Type{Array{T,n}})=Array{real(T),n}

eps(x...)=Base.eps(x...)
eps{T<:Real}(::Type{Complex{T}})=eps(real(T))
eps{T<:Real}(z::Complex{T})=eps(abs(z))
eps{T<:Number}(::Type{Vector{T}})=eps(T)
eps{k,T<:Number}(::Type{Vec{k,T}})=eps(T)


# BLAS

dotu(f::Vector{Complex{Float64}},g::Vector{Complex{Float64}})=BLAS.dotu(f,g)
dotu{N<:Real}(f::Vector{Complex{Float64}},g::Vector{N})=dot(conj(f),g)
dotu{N<:Real,T<:Number}(f::Vector{N},g::Vector{T})=dot(f,g)

# implement muladd default
muladd(a,b,c)=a*b+c
muladd(a::Number,b::Number,c::Number)=Base.muladd(a,b,c)


for TYP in (:Float64,:Float32,:Complex128,:Complex64)
    @eval scal!{T<:$TYP}(n::Integer,cst::$TYP,ret::DenseArray{T},k::Integer) =
            BLAS.scal!(n,cst,ret,k)
end

typealias BlasNumber Union{Float64,Float32,Complex128,Complex64}
scal!{T<:BlasNumber}(n::Integer,cst::BlasNumber,ret::DenseArray{T},k::Integer) =
    BLAS.scal!(n,T(cst),ret,k)

function scal!(n::Integer,cst::Number,ret::AbstractArray,k::Integer)
    @assert k*n ≤ length(ret)
    @simd for j=1:k:k*(n-1)+1
        @inbounds ret[j] *= cst
    end
    ret
end

scal!(cst::Number,v::AbstractArray) = scal!(length(v),cst,v,1)



## Helper routines
alternatingvector(n::Integer) = 2*mod([1:n],2) .- 1

function alternatesign!(v::Vector)
    n=length(v)
    for k=2:2:n
        v[k]=-v[k]
    end

    v
end

alternatesign(v::Vector)=alternatesign!(copy(v))

function alternatingsum(v::Vector)
    ret=zero(eltype(v))
    s=1
    @inbounds for k=1:length(v)
        ret+=s*v[k]
        s*=-1
    end

    ret
end



function pad!{T}(f::Vector{T},n::Integer)
	if n > length(f)
		append!(f,zeros(T,n - length(f)))
	else
		resize!(f,n)
	end
end


function pad{T}(f::Vector{T},n::Integer)
	if n > length(f)
	   ret=Array(T,n)
	   ret[1:length(f)]=f
	   for j=length(f)+1:n
	       ret[j]=zero(T)
	   end
       ret
	else
        f[1:n]
	end
end

function pad(f::Vector{Any},n::Integer)
	if n > length(f)
        Any[f...,zeros(n - length(f))...]
	else
        f[1:n]
	end
end

function pad(v::Vector,n::Integer,m::Integer)
    @assert m==1
    pad(v,n)
end

function pad(A::Matrix,n::Integer,m::Integer)
    T=eltype(A)
	if n <= size(A,1) && m <= size(A,2)
        A[1:n,1:m]
	elseif n==0 || m==0
	   Array(T,n,m)  #fixes weird julia bug when T==None
    else
        ret = Array(T,n,m)
        minn=min(n,size(A,1))
        minm=min(m,size(A,2))
        for k=1:minn,j=1:minm
            @inbounds ret[k,j]=A[k,j]
        end
        for k=minn+1:n,j=1:minm
            @inbounds ret[k,j]=zero(T)
        end
        for k=1:n,j=minm+1:m
            @inbounds ret[k,j]=zero(T)
        end
        for k=minn+1:n,j=minm+1:m
            @inbounds ret[k,j]=zero(T)
        end

        ret
	end
end

pad(A::Matrix,::Colon,m::Integer)=pad(A,size(A,1),m)
pad(A::Matrix,n::Integer,::Colon)=pad(A,n,size(A,2))


#TODO:padleft!

function padleft(f::Vector,n::Integer)
	if (n > length(f))
        [zeros(n - length(f)),f]
	else
        f[end-n+1:end]
	end
end



##chop!
function chop!(c::Vector,tol::Real)
    @assert tol >= 0

    for k=length(c):-1:1
        if abs(c[k]) > tol
            resize!(c,k)
            return c
        end
    end

    resize!(c,0)
    c
end

chop(f::Vector,tol)=chop!(copy(f),tol)
chop!(f::Vector)=chop!(f,eps())


function chop!(A::Array,tol)
    for k=size(A,1):-1:1
        if norm(A[k,:])>tol
            A=A[1:k,:]
            break
        end
    end
    for k=size(A,2):-1:1
        if norm(A[:,k])>tol
            A=A[:,1:k]
            break
        end
    end
    return A
end
chop(A::Array,tol)=chop!(A,tol)#replace by chop!(copy(A),tol) when chop! is actually in-place.



## interlace



function interlace(v::Union{Vector{Any},Tuple})
    #determine type
    T=Float64
    for vk in v
        if isa(vk,Vector{Complex{Float64}})
            T=Complex{Float64}
        end
    end
    b=Array(Vector{T},length(v))
    for k=1:length(v)
        b[k]=v[k]
    end
    interlace(b)
end

function interlace{S<:Number,V<:Number}(a::Vector{S},b::Vector{V})
    na=length(a);nb=length(b)
    T=promote_type(S,V)
    if nb≥na
        ret=zeros(T,2nb)
        ret[1:2:1+2*(na-1)]=a
        ret[2:2:end]=b
        ret
    else
        ret=zeros(T,2na-1)
        ret[1:2:end]=a
        if !isempty(b)
            ret[2:2:2+2*(nb-1)]=b
        end
        ret
    end
end

function interlace(a::Vector,b::Vector)
    na=length(a);nb=length(b)
    T=promote_type(eltype(a),eltype(b))
    if nb≥na
        ret=Array(T,2nb)
        ret[1:2:1+2*(na-1)]=a
        ret[2:2:end]=b
        ret
    else
        ret=Array(T,2na-1)
        ret[1:2:end]=a
        if !isempty(b)
            ret[2:2:2+2*(nb-1)]=b
        end
        ret
    end
end


# this limits the dimension of the padding to da and db
# interlacing every other coefficient until then

interlace(a::Union{Tuple,AbstractVector};dimensions=fill(Inf,length(a)))=dim_interlace(a,dimensions)
dim_interlace(a,dimensions)=
        dim_interlace(mapreduce(eltype,promote_type,a),a,dimensions)
function dim_interlace{T}(::Type{T},a,d)
    @assert length(d)==length(a)
    m=length(a)
    for j=1:m
        @assert length(a[j])≤d[j]
    end
    ret=Array(T,0)
    n=mapreduce(length,max,a)   # the max length
    for k=1:n, j=1:m
        if k ≤ length(a[j])
            push!(ret,a[j][k])
        elseif k ≤ d[j]
            # only add zero if we are less than the dimension dictated by d
            push!(ret,zero(T))
        end
    end

    ret
end





## svfft

##FFT That interlaces coefficients

plan_svfft(x::Vector) = plan_fft(x)
plan_isvfft(x::Vector) = plan_ifft(x)

function svfft(v::Vector,plan)
    n=length(v)
    v=plan*v/n
    if mod(n,2) == 0
        ind=div(n,2)
        v=alternatesign!(v)
        interlace(v[1:ind],
                  flipdim(v[ind+1:end],1))
    elseif mod(n,4)==3
        ind=div(n+1,2)
        interlace(alternatesign!(v[1:ind]),
                  -flipdim(alternatesign!(v[ind+1:end]),1))
    else #mod(length(v),4)==1
        ind=div(n+1,2)
        interlace(alternatesign!(v[1:ind]),
                  flipdim(alternatesign!(v[ind+1:end]),1))
    end
end

function isvfft(sv::Vector,plan)
    n=length(sv)

    if mod(n,2) == 0
        v=alternatesign!([sv[1:2:end];flipdim(sv[2:2:end],1)])
    elseif mod(n,4)==3
        v=[alternatesign!(sv[1:2:end]);
           -alternatesign!(flipdim(sv[2:2:end],1))]
    else #mod(length(v),4)==1
        v=[alternatesign!(sv[1:2:end]);
           alternatesign!(flipdim(sv[2:2:end],1))]
    end

    plan*(n*v)
end




## slnorm gives the norm of a slice of a matrix

function slnorm{T}(u::AbstractArray{T},r::Range,::Colon)
    ret = zero(real(T))
    for k=r
        @simd for j=1:size(u,2)
            #@inbounds
            ret=max(norm(u[k,j]),ret)
        end
    end
    ret
end


function slnorm(m::AbstractMatrix,kr::Range,jr::Range)
    ret=0.0
    for j=jr
        for k=kr
            @inbounds ret=ret+abs2(m[k,j])
        end
    end
    ret
end

slnorm(m::Matrix,kr::Range,jr::Integer)=slnorm(m,kr,jr:jr)
slnorm(m::Matrix,kr::Integer,jr::Range)=slnorm(m,kr:kr,jr)


function slnorm{T}(B::BandedMatrix{T},r::Range,::Colon)
    ret = zero(real(T))
    m=size(B,2)
    for k=r
        @simd for j=max(1,k-B.l):min(k+B.u,m)
            #@inbounds
            ret=max(norm(B[k,j]),ret)
        end
    end
    ret
end






## New Inf

const ∞ = Irrational{:∞}()

Base.show(io::IO, x::Irrational{:∞}) = print(io, "∞")
Base.convert{F<:AbstractFloat}(::Type{F},::Irrational{:∞}) = convert(F,Inf)

## My Count

abstract AbstractCount{S<:Number}

immutable UnitCount{S<:Number} <: AbstractCount{S}
    start::S
end

immutable Count{S<:Number} <: AbstractCount{S}
    start::S
    step::S
end
countfrom(start::Number, step::Number) = Count(promote(start, step)...)
countfrom(start::Number)               = UnitCount(start)
countfrom()                            = UnitCount(1)


Base.eltype{S}(::Type{AbstractCount{S}}) = S
Base.eltype{AS<:AbstractCount}(::Type{AS}) = eltype(super(AS))

Base.step(it::Count) = it.step
Base.step(it::UnitCount) = 1

Base.start(it::AbstractCount) = it.start
Base.next(it::AbstractCount, state) = (state, state + step(it))
Base.done(it::AbstractCount, state) = false

Base.length(it::AbstractCount) = ∞

getindex(it::Count,k) = it.start + it.step*(k-1)
getindex(it::UnitCount,k) = it.start + k - 1


Base.colon(a::Real,b::Irrational{:∞}) = countfrom(a)
Base.colon(::Irrational{:∞},::AbstractFloat,::Irrational{:∞}) = [∞]
Base.colon(a::Real,st::Real,b::Irrational{:∞}) = countfrom(a,st)

Base.isinf(::Irrational{:∞}) = true
