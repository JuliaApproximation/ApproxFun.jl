export dotu
import Base: chop

# Used for spaces not defined yet
struct UnsetNumber <: Number  end
promote_rule(::Type{UnsetNumber},::Type{N}) where {N<:Number} = N

# Test the number of arguments a function takes
hasnumargs(f,k) = applicable(f,zeros(k)...)


# fast implementation of isapprox with atol a non-keyword argument in most cases
isapprox_atol(a,b,atol;kwds...) = isapprox(a,b;atol=atol,kwds...)
isapprox_atol(a::Vec,b::Vec,atol::Real=0;kwds...) = isapprox_atol(collect(a),collect(b),atol;kwds...)
function isapprox_atol(x::Number, y::Number, atol::Real=0; rtol::Real=Base.rtoldefault(x,y))
    x == y || (isfinite(x) && isfinite(y) && abs(x-y) <= atol + rtol*max(abs(x), abs(y)))
end
function isapprox_atol(x::AbstractArray{T}, y::AbstractArray{S},atol::Real=0; rtol::Real=Base.rtoldefault(T,S), norm::Function=vecnorm) where {T<:Number,S<:Number}
    d = norm(x - y)
    if isfinite(d)
        return d <= atol + rtol*max(norm(x), norm(y))
    else
        # Fall back to a component-wise approximate comparison
        return all(ab -> isapprox(ab[1], ab[2]; rtol=rtol, atol=atol), zip(x, y))
    end
end

# The second case handles zero
isapproxinteger(x) = isapprox(x,round(Int,x))  || isapprox(x+1,round(Int,x+1))


# This creates ApproxFun.real, ApproxFun.eps and ApproxFun.dou
# which we override for default julia types
real(x...) = Base.real(x...)
real(::Type{UnsetNumber}) = UnsetNumber
real(::Type{Array{T,n}}) where {T<:Real,n} = Array{T,n}
real(::Type{Array{T,n}}) where {T<:Complex,n} = Array{real(T),n}
real(::Type{Vec{N,T}}) where {N,T<:Real} = Vec{N,T}
real(::Type{Vec{N,T}}) where {N,T<:Complex} = Vec{N,real(T)}


eps(x...) = Base.eps(x...)
eps(x) = Base.eps(x)

eps(::Type{T}) where {T<:Integer} = zero(T)

eps(::Type{Complex{T}}) where {T<:Real} = eps(real(T))
eps(z::Complex{T}) where {T<:Real} = eps(abs(z))
eps(::Type{Dual{Complex{T}}}) where {T<:Real} = eps(real(T))
eps(z::Dual{Complex{T}}) where {T<:Real} = eps(abs(z))

eps(::Type{Vector{T}}) where {T<:Number} = eps(T)
eps(::Type{Vec{k,T}}) where {k,T<:Number} = eps(T)


isnan(x) = Base.isnan(x)
isnan(x::Vec) = map(isnan,x)


# BLAS


# implement muladd default
muladd(a,b,c) = a*b+c
muladd(a::Number,b::Number,c::Number) = Base.muladd(a,b,c)


for TYP in (:Float64,:Float32,:ComplexF64,:ComplexF32)
    @eval scal!(n::Integer,cst::$TYP,ret::DenseArray{T},k::Integer) where {T<:$TYP} =
            BLAS.scal!(n,cst,ret,k)
end


scal!(n::Integer,cst::BlasFloat,ret::DenseArray{T},k::Integer) where {T<:BlasFloat} =
    BLAS.scal!(n,T(cst),ret,k)

function scal!(n::Integer,cst::Number,ret::AbstractArray,k::Integer)
    @assert k*n ≤ length(ret)
    @simd for j=1:k:k*(n-1)+1
        @inbounds ret[j] *= cst
    end
    ret
end

scal!(cst::Number,v::AbstractArray) = scal!(length(v),cst,v,1)



# Helper routines

function reverseeven!(x::AbstractVector)
    n = length(x)
    if iseven(n)
        @inbounds @simd for k=2:2:n÷2
            x[k],x[n+2-k] = x[n+2-k],x[k]
        end
    else
        @inbounds @simd for k=2:2:n÷2
            x[k],x[n+1-k] = x[n+1-k],x[k]
        end
    end
    x
end

function negateeven!(x::AbstractVector)
    @inbounds @simd for k = 2:2:length(x)
        x[k] *= -1
    end
    x
end

#checkerboard, same as applying negativeeven! to all rows then all columns
function negateeven!(X::AbstractMatrix)
    for j = 1:2:size(X,2)
        @inbounds @simd for k = 2:2:size(X,1)
            X[k,j] *= -1
        end
    end
    for j = 2:2:size(X,2)
        @inbounds @simd for k = 1:2:size(X,1)
            X[k,j] *= -1
        end
    end
    X
end

const alternatesign! = negateeven!

alternatesign(v::AbstractVector) = alternatesign!(copy(v))

alternatingvector(n::Integer) = 2*mod([1:n],2) .- 1

function alternatingsum(v::AbstractVector)
    ret = zero(eltype(v))
    s = 1
    @inbounds for k=1:length(v)
        ret+=s*v[k]
        s*=-1
    end

    ret
end

# Sum Hadamard product of vectors up to minimum over lengths
function mindotu(a::AbstractVector,b::AbstractVector)
    ret,m = zero(promote_type(eltype(a),eltype(b))),min(length(a),length(b))
    @inbounds @simd for i=m:-1:1 ret += a[i]*b[i] end
    ret
end


# efficiently resize a Matrix.  Note it doesn't change the input ptr
function unsafe_resize!(W::AbstractMatrix,::Colon,m::Integer)
    if m == size(W,2)
        W
    else
        n=size(W,1)
        reshape(resize!(vec(W),n*m),n,m)
    end
end

function unsafe_resize!(W::AbstractMatrix,n::Integer,::Colon)
    N=size(W,1)
    if n == N
        W
    elseif n < N
        W[1:n,:]
    else
        m=size(W,2)
        ret=Matrix{eltype(W)}(n,m)
        ret[1:N,:] = W
        ret
    end
end

function unsafe_resize!(W::AbstractMatrix,n::Integer,m::Integer)
    N=size(W,1)
    if n == N
        unsafe_resize!(W,:,m)
    else
        unsafe_resize!(unsafe_resize!(W,n,:),:,m)
    end
end


function pad!(f::AbstractVector{T},n::Integer) where T
	if n > length(f)
		append!(f,zeros(T,n - length(f)))
	else
		resize!(f,n)
	end
end


function pad(f::AbstractVector{T},n::Integer) where T
	if n > length(f)
	   ret=Vector{T}(undef, n)
	   ret[1:length(f)]=f
	   for j=length(f)+1:n
	       ret[j]=zero(T)
	   end
       ret
	else
        f[1:n]
	end
end

function pad(f::AbstractVector{Any},n::Integer)
	if n > length(f)
        Any[f...,zeros(n - length(f))...]
	else
        f[1:n]
	end
end

function pad(v::AbstractVector,n::Integer,m::Integer)
    @assert m==1
    pad(v,n)
end

function pad(A::AbstractMatrix,n::Integer,m::Integer)
    T=eltype(A)
	if n <= size(A,1) && m <= size(A,2)
        A[1:n,1:m]
	elseif n==0 || m==0
	   Matrix{T}(undef,n,m)  #fixes weird julia bug when T==None
    else
        ret = Matrix{T}(undef,n,m)
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

pad(A::AbstractMatrix,::Colon,m::Integer) = pad(A,size(A,1),m)
pad(A::AbstractMatrix,n::Integer,::Colon) = pad(A,n,size(A,2))



#TODO:padleft!

function padleft(f::AbstractVector,n::Integer)
	if (n > length(f))
        [zeros(n - length(f)),f]
	else
        f[end-n+1:end]
	end
end



##chop!
function chop!(c::AbstractVector,tol::Real)
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

chop(f::AbstractVector,tol) = chop!(copy(f),tol)
chop!(f::AbstractVector) = chop!(f,eps())


function chop!(A::AbstractArray,tol)
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
chop(A::AbstractArray,tol)=chop!(A,tol)#replace by chop!(copy(A),tol) when chop! is actually in-place.



## interlace



function interlace(v::Union{Vector{Any},Tuple})
    #determine type
    T=Float64
    for vk in v
        if isa(vk,Vector{Complex{Float64}})
            T=Complex{Float64}
        end
    end
    b=Vector{Vector{T}}(length(v))
    for k=1:length(v)
        b[k]=v[k]
    end
    interlace(b)
end

function interlace(a::AbstractVector{S},b::AbstractVector{V}) where {S<:Number,V<:Number}
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

function interlace(a::AbstractVector,b::AbstractVector)
    na=length(a);nb=length(b)
    T=promote_type(eltype(a),eltype(b))
    if nb≥na
        ret=Vector{T}(2nb)
        ret[1:2:1+2*(na-1)]=a
        ret[2:2:end]=b
        ret
    else
        ret=Vector{T}(2na-1)
        ret[1:2:end]=a
        if !isempty(b)
            ret[2:2:2+2*(nb-1)]=b
        end
        ret
    end
end


### In-place O(n) interlacing

function highestleader(n::Int)
    i = 1
    while 3i < n i *= 3 end
    i
end

function nextindex(i::Int,n::Int)
    i <<= 1
    while i > n
        i -= n + 1
    end
    i
end

function cycle_rotate!(v::AbstractVector, leader::Int, it::Int, twom::Int)
    i = nextindex(leader, twom)
    while i != leader
        idx1, idx2 = it + i - 1, it + leader - 1
        @inbounds v[idx1], v[idx2] = v[idx2], v[idx1]
        i = nextindex(i, twom)
    end
    v
end

function right_cyclic_shift!(v::AbstractVector, it::Int, m::Int, n::Int)
    itpm = it + m
    itpmm1 = itpm - 1
    itpmpnm1 = itpmm1 + n
    reverse!(v, itpm, itpmpnm1)
    reverse!(v, itpm, itpmm1 + m)
    reverse!(v, itpm + m, itpmpnm1)
    v
end

"""
This function implements the algorithm described in:

    P. Jain, "A simple in-place algorithm for in-shuffle," arXiv:0805.1598, 2008.
"""
function interlace!(v::AbstractVector,offset::Int)
    N = length(v)
    if N < 2 + offset
        return v
    end

    it = 1 + offset
    m = 0
    n = 1

    while m < n
        twom = N + 1 - it
        h = highestleader(twom)
        m = h > 1 ? h÷2 : 1
        n = twom÷2

        right_cyclic_shift!(v,it,m,n)

        leader = 1
        while leader < 2m
            cycle_rotate!(v, leader, it, 2m)
            leader *= 3
        end

        it += 2m
    end
    v
end

## slnorm gives the norm of a slice of a matrix

function slnorm(u::AbstractMatrix,r::AbstractRange,::Colon)
    ret = zero(real(eltype(u)))
    for k=r
        @simd for j=1:size(u,2)
            #@inbounds
            ret=max(norm(u[k,j]),ret)
        end
    end
    ret
end


function slnorm(m::AbstractMatrix,kr::AbstractRange,jr::AbstractRange)
    ret=zero(real(eltype(m)))
    for j=jr
        nrm=zero(real(eltype(m)))
        for k=kr
            @inbounds nrm+=abs2(m[k,j])
        end
        ret=max(sqrt(nrm),ret)
    end
    ret
end

slnorm(m::AbstractMatrix,kr::AbstractRange,jr::Integer) = slnorm(m,kr,jr:jr)
slnorm(m::AbstractMatrix,kr::Integer,jr::AbstractRange) = slnorm(m,kr:kr,jr)


function slnorm(B::BandedMatrix{T},r::AbstractRange,::Colon) where T
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


slnorm(m::AbstractMatrix,k::Integer,::Colon) = slnorm(m,k,1:size(m,2))
slnorm(m::AbstractMatrix,::Colon,j::Integer) = slnorm(m,1:size(m,1),j)



## New Inf

# angle is π*a where a is (false==0) and (true==1)
struct Infinity{T}
    angle::T
end

Infinity() = Infinity(false)
const ∞ = Infinity()


isinf(::Infinity) = true
isfinite(::Infinity) = false
sign(y::Infinity{B}) where {B<:Integer} = mod(y.angle,2)==0 ? 1 : -1
angle(x::Infinity) = π*x.angle

function show(io::IO, y::Infinity{B}) where B<:Integer
    if sign(y) == 1
        print(io, "∞")
    else
        print(io, "-∞")
    end
end

show(io::IO,x::Infinity) = print(io,"$(exp(im*π*x.angle))∞")

==(x::Infinity,y::Infinity) = x.angle == y.angle
for TYP in (:Dual,:Number)
    @eval begin
        ==(x::Infinity,y::$TYP) = isinf(y) && angle(y) == angle(x)
        ==(y::$TYP,x::Infinity) = x == y
    end
end
isless(x::Infinity{Bool}, y::Infinity{Bool}) = x.angle && !y.angle
isless(x::Number, y::Infinity{Bool}) = !y.angle && x ≠ ∞
isless(x::Infinity{Bool}, y::Number) = x.angle && y ≠ -∞
isless(x::Block{1}, y::Infinity{Bool}) = isless(Int(x), y)
isless(x::Infinity{Bool}, y::Block{1}) = isless(x, Int(y))

-(y::Infinity{B}) where {B<:Integer} = sign(y)==1 ? Infinity(one(B)) : Infinity(zero(B))

function +(x::Infinity{B}, y::Infinity{B}) where B
    if x.angle != y.angle
        error("Angles must be the same to add ∞")
    end
    x
end

for T in (:BlasFloat,:Integer,:(Complex{Int}))
    @eval begin
        +(::$T,y::Infinity) = y
        +(y::Infinity,::$T) = y
        -(y::Infinity,::$T) = y
        -(::$T,y::Infinity) = -y
    end
end


# ⊻ is xor
*(a::Infinity{Bool},b::Infinity{Bool}) = Infinity(a.angle ⊻ b.angle)
*(a::Infinity,b::Infinity) = Infinity(a.angle + b.angle)

for T in (:Dual,:Bool,:Integer,:AbstractFloat)
    @eval begin
        *(a::$T,y::Infinity) = a>0 ? y : (-y)
        *(y::Infinity,a::$T) = a*y
    end
end

*(a::Number,y::Infinity) = Infinity(y.angle+angle(a)/π)
*(y::Infinity,a::Number) = a*y

for OP in (:fld,:cld,:div)
  @eval $OP(y::Infinity,a::Number) = y*(1/sign(a))
end

min(x::Infinity{B},y::Infinity{B}) where {B<:Integer} = sign(x)==-1 ? x : y
max(x::Infinity{B},::Infinity{B}) where {B<:Integer} = sign(x)==1 ? x : y
min(x::Real,y::Infinity{B}) where {B<:Integer} = sign(y)==1 ? x : y
min(x::Infinity{B},y::Real) where {B<:Integer} = min(y,x)
max(x::Real,y::Infinity{B}) where {B<:Integer} = sign(y)==1 ? y : x
max(x::Infinity{B},y::Real) where {B<:Integer} = max(y,x)

for OP in (:<,:<=)
    @eval begin
        $OP(x::Real,y::Infinity{B}) where {B<:Integer} = sign(y)==1
        $OP(y::Infinity{B},x::Real) where {B<:Integer} = sign(y)==-1
    end
end

for OP in (:>,:>=)
    @eval begin
        $OP(x::Real,y::Infinity{B}) where {B<:Integer} = sign(y)==-1
        $OP(y::Infinity{B},x::Real) where {B<:Integer} = sign(y)==1
    end
end


abstract type Iterator end


function findfirst(testf::Function, A::Iterator)
    for (k,v) in enumerate(A)
        testf(v) && return k
    end
    return 0
end

# Take -- iterate through the first n elements

struct Take{I,T} <: AbstractVector{T}
    xs::I
    n::Int
    function Take{I,T}(xs,n) where {I,T}
        @assert n ≤ length(xs)
        new{I,T}(xs,n)
    end
end

Take(xs,n) = Take{typeof(xs),eltype(xs)}(xs,n)

"""
    take(iter, n)

An iterator that generates at most the first `n` elements of `iter`.

```jldoctest
julia> a = 1:2:11
1:2:11

julia> collect(a)
6-element Array{Int64,1}:
  1
  3
  5
  7
  9
 11

julia> collect(take(a,3))
3-element Array{Int64,1}:
 1
 3
 5
```
"""
take(xs, n::Int) = Take(xs, n)

eltype(::Type{Take{I}}) where {I} = eltype(I)
length(t::Take) = t.n
size(t::Take) = (length(t),)
start(it::Take) = (it.n, start(it.xs))
function next(it::Take, state)
    n, xs_state = state
    v, xs_state = next(it.xs, xs_state)
    return v, (n - 1, xs_state)
end

function done(it::Take, state)
    n, xs_state = state
    return n <= 0 || done(it.xs, xs_state)
end


function getindex(it::Take,k)
    !isempty(k) && maximum(k) > it.n && throw(BoundsError())

    it.xs[k]
end

getindex(it::Take,k::CartesianIndex{1}) = it[k[1]]

function sum(it::Take)
    ret = zero(eltype(it))
    for a in it
        ret += a
    end
    ret
end

pad(it::Take,n::Integer) = pad!(collect(it),n)



# Re-implementation of Base iterators
# to use ∞ and allow getindex

abstract type AbstractRepeated{T} <: Iterator end

eltype(::Type{AbstractRepeated{T}}) where {T} = T
eltype(::Type{R}) where {R<:AbstractRepeated} = eltype(super(R))
eltype(::AbstractRepeated{T}) where {T} = T

step(::AbstractRepeated) = 0

start(::AbstractRepeated) = nothing
next(it::AbstractRepeated,state) = value(it),nothing
done(::AbstractRepeated,state) = false

length(::AbstractRepeated) = ∞

getindex(it::AbstractRepeated,k::Integer) = value(it)
getindex(it::AbstractRepeated,k::AbstractRange) = take(it,length(k))


maximum(r::AbstractRepeated) = value(r)
minimum(r::AbstractRepeated) = value(r)

sum(it::Take{AR}) where {AR<:AbstractRepeated} = it.n*value(it.xs)

struct ZeroRepeated{T} <: AbstractRepeated{T} end

ZeroRepeated(::Type{T}) where {T} = ZeroRepeated{T}()

value(::ZeroRepeated{T}) where {T} = zero(T)
sum(r::ZeroRepeated) = value(r)

struct Repeated{T} <: AbstractRepeated{T}
    x::T
    function Repeated{T}(x::T) where T
        # TODO: Add ZeroRepeated type.
        if x == zero(T)
            error("Zero repeated not supported to maintain type stability")
        end

        new{T}(x)
    end
end

Repeated(x) = Repeated{typeof(x)}(x)


value(r::Repeated) = r.x

sum(r::Repeated) = r.x > 0 ? ∞ : -∞



function repeated(x)
    if x == zero(x)
        error("Use ZeroRepeated to repeat zeros")
    end
    Repeated(x)
end
repeated(x,::Infinity{Bool}) = repeated(x)
repeated(x,m::Integer) = take(repeated(x),m)


abstract type AbstractCount{S} <: Iterator end

struct UnitCount{S} <: AbstractCount{S}
    start::S
end

struct Count{S} <: AbstractCount{S}
    start::S
    step::S
end


countfrom(start::Number, step::Number) = Count(promote(start, step)...)
countfrom(start::Number)               = UnitCount(start)
countfrom()                            = UnitCount(1)


eltype(::Type{AbstractCount{S}}) where {S} = S
eltype(::Type{AS}) where {AS<:AbstractCount} = eltype(supertype(AS))
eltype(::AbstractCount{S}) where {S} = S

step(it::Count) = it.step
step(it::UnitCount) = 1

start(it::AbstractCount) = it.start
next(it::AbstractCount, state) = (state, state + step(it))
done(it::AbstractCount, state) = false

length(it::AbstractCount) = ∞

getindex(it::Count,k) = it.start + it.step*(k-1)
getindex(it::UnitCount,k) = (it.start-1) + k
getindex(it::AbstractRepeated,k::AbstractCount) = it

# use reindex, copied from Base
reindex(V, idxs::Tuple{AbstractCount, Vararg{Any}}, subidxs::Tuple{Any, Vararg{Any}}) =
    (Base.@_propagate_inbounds_meta; (idxs[1][subidxs[1]], reindex(V, tail(idxs), tail(subidxs))...))

# special functions
maximum(it::UnitCount{S}) where {S<:Real} = ∞
minimum(it::UnitCount{S}) where {S<:Real} = it.start

maximum(it::Count{S}) where {S<:Real} = it.step > 0 ? ∞ : it.start
minimum(it::Count{S}) where {S<:Real} = it.step < 0 ? -∞ : it.start

sum(it::UnitCount) = ∞
sum(it::Count) = it.step > 0 ? ∞ : -∞

@inline first(it::AbstractCount) = it.start

last(it::UnitCount{S}) where {S<:Real} = ∞
last(it::Count{S}) where {S<:Real} =
    it.step > 0 ? ∞ : (it.step < 0 ? -∞ : error("zero step not supported"))


function (:)(a::Real,b::Infinity{Bool})
    if b.angle
        throw(ArgumentError("Cannot create $a:-∞"))
    end
    countfrom(a)
end
function (:)(a::Infinity{Bool},st::AbstractFloat,b::Infinity{Bool})
    if a ≠ b
        throw(ArgumentError("Cannot create $a:$st:$b"))
    end
    [a]
end
function (:)(a::Real,st::Real,b::Infinity{Bool})
    if st == 0
        throw(ArgumentError("step cannot be zero"))
    elseif b.angle == st > 0
        throw(ArgumentError("Cannot create $a:$st:$b"))
    else
        countfrom(a,st)
    end
end



intersect(a::UnitCount, b::UnitCount) = UnitCount(max(first(a), first(b)))
intersect(a::AbstractCount, b::AbstractCount) = error("Not implemented")

intersect(a::UnitCount, b::AbstractRange) = intersect(first(a):last(b), b)
intersect(a::AbstractRange, b::UnitCount) = intersect(a, first(b):last(a))

intersect(a::Count, b::AbstractRange) = intersect(first(a):step(a):last(b), b)
intersect(a::AbstractRange, b::Count) = intersect(a, first(b):step(b):last(a))



struct CumSumIterator{CC} <: Iterator
    iterator::CC
end

eltype(::Type{CumSumIterator{S}}) where {S} = eltype(S)
eltype(CC::CumSumIterator) = eltype(CC.iterator)

start(it::CumSumIterator) = (0,start(it.iterator))
function next(it::CumSumIterator, state)
    a,nx_st=next(it.iterator,state[2])
    cs=state[1]+a
    (cs,(cs,nx_st))
end
done(it::CumSumIterator, state) = done(it.iterator,state[2])

length(it::CumSumIterator) = length(it.iterator)

getindex(it::CumSumIterator{AC},k) where {AC<:UnitCount} = it.iterator.start*k + ((k*(k-1))÷2)
getindex(it::CumSumIterator{AC},k) where {AC<:Count} = it.iterator.start*k + step(it.iterator)*((k*(k-1))÷2)





## BandedMatrix



pad!(A::BandedMatrix,n,::Colon) = pad!(A,n,n+A.u)  # Default is to get all columns
columnrange(A,row::Integer) = max(1,row+bandinds(A,1)):row+bandinds(A,2)



## Store iterator
mutable struct CachedIterator{T,IT,ST} <: Iterator
    iterator::IT
    storage::Vector{T}
    state::ST
    length::Int

    CachedIterator{T,IT,ST}(it::IT) where {T,IT,ST} = new{T,IT,ST}(it,Vector{T}(),start(it),0)
end

CachedIterator(it) = CachedIterator{eltype(it),typeof(it),typeof(start(it))}(it)

function resize!(it::CachedIterator,n::Integer)
    m = it.length
    if n > m
        if n > length(it.storage)
            resize!(it.storage,2n)
        end

        for k = m+1:n
            if done(it.iterator,it.state)
                it.length = k-1
                return it
            end
            val,it.state = next(it.iterator,it.state)
            @inbounds it.storage[k] = val
        end

        it.length = n
    end
    it
end


eltype(it::CachedIterator{T}) where {T} = T
start(it::CachedIterator) = 1
next(it::CachedIterator,st::Int) = (it[st],st+1)
done(it::CachedIterator,st::Int) = st == it.length + 1 &&
                                        done(it.iterator,it.state)

function getindex(it::CachedIterator,k)
    mx = maximum(k)
    if mx > length(it) || mx < 1
        throw(BoundsError(it,k))
    end
    resize!(it,isempty(k) ? 0 : mx).storage[k]
end
function findfirst(f::Function,A::CachedIterator)
    k=1
    for c in A
        if f(c)
            return k
        end
        k+=1
    end
    return 0
end

function findfirst(A::CachedIterator,x)
    k=1
    for c in A
        if c == x
            return k
        end
        k+=1
    end
    return 0
end

length(A::CachedIterator) = length(A.iterator)


# The following don't need caching
cache(A::AbstractVector{T}) where {T<:Number} = A
cache(A::AbstractRange) = A
cache(A::AbstractCount) = A



## From Julia v0.5 code


# flatten an iterator of iterators
# we add indexing

struct Flatten{I}
    it::I
end

"""
    flatten(iter)

Given an iterator that yields iterators, return an iterator that yields the
elements of those iterators.
Put differently, the elements of the argument iterator are concatenated. Example:

    julia> collect(flatten((1:2, 8:9)))
    4-element Array{Int64,1}:
     1
     2
     8
     9
"""
flatten(itr) = Flatten(itr)

eltype(f::Flatten) = mapreduce(eltype,promote_type,f.it)


start(f::Flatten) = 1, map(start,f.it)

function next(f::Flatten, state)
    k, sts = state
    if !done(f.it[k],sts[k])
        a,nst = next(f.it[k],sts[k])
        a, (k,(sts[1:k-1]...,nst,sts[k+1:end]...))
    else
        next(f,(k+1,sts))
    end
end


@inline done(f::Flatten, state) =
    state[1] == length(f) && done(f.it[end],state[2][end])

length(f::Flatten) = mapreduce(length,+,f.it)

eachindex(f::Flatten) = 1:length(f)

function getindex(f::Flatten,k::Int)
    sh = 0
    for it in f.it
        n = length(it)
        if k ≤ sh + n
            return it[k-sh]
        else
            sh += n
        end
    end

    throw(BoundsError())
end

getindex(f::Flatten,kr::UnitRange{Int}) = eltype(f)[f[k] for k in kr]

sum(f::Flatten) = mapreduce(sum,+,f.it)


maximum(f::Flatten) = mapreduce(maximum,max,f.it)
minimum(f::Flatten) = mapreduce(minimum,min,f.it)


## Iterator Algebra

broadcast(op,f::Flatten,c...) = Flatten(map(it->op(it,c...),f.it))


broadcast(op,a::AbstractRepeated,b::AbstractRepeated) = repeated(op.(value(a),value(b)))
broadcast(op,a::AbstractRepeated,b::Number) = repeated(op.(value(a),b))
broadcast(op,a::Number,b::AbstractRepeated) = repeated(op.(a,value(b)))

broadcast(op,a::AbstractCount,b::AbstractRepeated) = op.(a,value(b))
broadcast(op,a::AbstractRepeated,b::AbstractCount) = op.(value(a),b)

function broadcast(op,a::Flatten,b::AbstractRepeated)
    @assert isinf(length(a.it[end]))
    flatten(map(it->op.(it,value(b)),a.it))
end
function broadcast(op,a::AbstractRepeated,b::Flatten)
    @assert isinf(length(b.it[end]))
    flatten(map(it->op.(value(a),it),b.it))
end
function broadcast(op,a::Flatten,b::AbstractCount)
    K=0
    it=tuple()
    for k=1:length(a.it)
        it=(it...,op.(a.it[k],b[K+1:K+length(a.it[k])]))
        K+=length(a.it[k])
    end
    flatten(it)
end
function broadcast(op,a::AbstractCount,b::Flatten)
    K=0
    it=tuple()
    for k=1:length(b.it)
        it=(it...,op.(a[K+1:K+length(it[k])],b.it[k]))
        K+=length(b.it[k])
    end
    flatten(it)
end
function broadcast(op,a::Take,b::Take)
    n = length(a)
    @assert n == length(b)
    take(op.(a.xs,b.xs),n)
end
function broadcast(op,a::Take,b::Number)
    n = length(a)
    take(op.(a.xs,b),n)
end
function broadcast(op,a::Number,b::Take)
    n = length(b)
    take(op.(a,b.xs),n)
end

const InfiniteIterators = Union{AbstractRepeated,AbstractCount,Flatten}

+(a::InfiniteIterators) = a
-(a::ZeroRepeated) = a
-(a::Repeated) = Repeated(-value(a))
-(a::Flatten) = Flatten(map(-,a.it))

for OP in (:+,:-)
    @eval begin
        $OP(a::AbstractRepeated,b::AbstractRepeated) = repeated($OP(value(a),value(b)))
        $OP(a::Number,b::AbstractRepeated) = repeated($OP(a,value(b)))
        $OP(a::AbstractRepeated,b::Number) = repeated($OP(value(a),b))

        $OP(a::AbstractCount,b::AbstractRepeated) = $OP(a,value(b))
        $OP(a::AbstractRepeated,b::AbstractCount) = $OP(value(a),b)

        function $OP(a::Flatten,b::AbstractRepeated)
            @assert isinf(length(a.it[end]))
            flatten(map(it->$OP(it,value(b)),a.it))
        end
        function $OP(a::AbstractRepeated,b::Flatten)
            @assert isinf(length(b.it[end]))
            flatten(map(it->$OP(value(a),it),b.it))
        end

        function $OP(a::Flatten,b::AbstractCount)
            K=0
            it=tuple()
            for k=1:length(a.it)
                it=(it...,$OP(a.it[k],b[K+1:K+length(a.it[k])]))
                K+=length(a.it[k])
            end
            flatten(it)
        end
        function $OP(a::AbstractCount,b::Flatten)
            K=0
            it=tuple()
            for k=1:length(b.it)
                it=(it...,$OP(a[K+1:K+length(it[k])],b.it[k]))
                K+=length(b.it[k])
            end
            flatten(it)
        end

        function $OP(a::Take,b::Take)
            n = length(a)
            @assert n == length(b)
            take($OP(a.xs,b.xs),n)
        end

        function $OP(a::Take,b::Bool)
            n = length(a)
            take($OP(a.xs,b),n)
        end
        function $OP(a::Bool,b::Take)
            n = length(b)
            take($OP(a,b.xs),n)
        end

        function $OP(a::Take,b::Number)
            n = length(a)
            take($OP(a.xs,b),n)
        end
        function $OP(a::Number,b::Take)
            n = length(b)
            take($OP(a,b.xs),n)
        end
    end
end


for OP in (:+,:-)
    @eval begin
        $OP(a::ZeroRepeated,b::ZeroRepeated) = a
        $OP(a::Number,b::ZeroRepeated) = repeated(a)
        $OP(a::ZeroRepeated,b::Number) = repeated($OP(b))


        $OP(a::Flatten,b::ZeroRepeated) = a
        $OP(a::ZeroRepeated,b::Flatten) = $OP(b)

        $OP(a::AbstractCount,b::AbstractCount) =
            Count($OP(start(a),start(b)),$OP(step(a),step(b)))
        $OP(a::UnitCount,b::Number) = UnitCount($OP(a.start,b))
        $OP(a::Count,b::Number) = Count($OP(a.start,b),a.step)

        $OP(a::Number,b::AbstractCount) = $OP(b) + a

        broadcast(::typeof($OP),a::InfiniteIterators,b::InfiniteIterators) = $OP(a,b)
        broadcast(::typeof($OP),a::Number,b::InfiniteIterators) = $OP(a,b)
        broadcast(::typeof($OP),a::InfiniteIterators,b::Number) = $OP(a,b)
    end
end

broadcast(::typeof(*),a::ZeroRepeated,b::ZeroRepeated) = a
broadcast(::typeof(*),a::Number,b::AbstractCount) = Count(start(b)*a,step(b)*a)
broadcast(::typeof(*),b::AbstractCount,a::Number) = a*b


+(a::Number,b::UnitCount) = UnitCount(a+b.start)
+(a::Number,b::Count) = Count(a+b.start,b.step)
-(a::Number,b::UnitCount) = Count(a-b.start,-1)
-(a::Number,b::Count) = Count(a-b.start,-b.step)

*(a::Number,b::AbstractCount) = Count(a*start(b),a*step(b))
*(a::AbstractCount,b::Number) = b*a

function +(a::Flatten,b::Flatten)
    if isempty(a)
        @assert isempty(b)
        a
    elseif length(a.it) == 1
        a.it[1]+b
    elseif length(b.it) == 1
        a+b.it[1]
    elseif length(a.it[1]) == length(b.it[1])
        flatten((a.it[1]+b.it[1],(flatten(a.it[2:end])+flatten(b.it[2:end])).it...))
    elseif length(a.it[1]) < length(b.it[1])
        n=length(a.it[1])
        flatten((a.it[1]+b.it[1][1:n],
            (flatten(a.it[2:end])+flatten((b.it[1][n+1:end],b.it[2:end]...))).it...))
    else #length(a.it[1]) > length(b.it[1])
        n=length(a.it[2])
        flatten((a.it[1][1:n]+b.it[1],
            (flatten((a.it[1][n+1:end],b.it[2:end]...))+flatten(b.it[2:end])).it...))
    end
end


cumsum(r::Repeated) = r.x:r.x:(r.x>0 ? ∞ : -∞)
cumsum(r::Repeated{Bool}) = 1:∞
cumsum(r::ZeroRepeated) = r
cumsum(r::Iterator) = CumSumIterator(r)





function cumsum(f::Flatten)
    cs=zero(eltype(f))
    its = Vector{eltype(f.it)}(undef, 0)
    for it in f.it[1:end-1]
        c = cumsum(cs .+ it)
        push!(its,c)
        cs=last(c)
    end

    c=cs+cumsum(f.it[end])
    push!(its,c)
    Flatten(tuple(its...))
end


function pad(v,::Infinity{Bool})
    if isinf(length(v))
        v
    else
        flatten((v,ZeroRepeated(Int)))
    end
end

function pad(v::AbstractArray,::Infinity{Bool})
    if isinf(length(v))
        v
    else
        flatten((v,ZeroRepeated(Int)))
    end
end


## nocat
vnocat(A...) = Base.vect(A...)
hnocat(A...) = Base.typed_hcat(mapreduce(typeof,promote_type,A),A...)
hvnocat(rows,A...) = Base.typed_hvcat(mapreduce(typeof,promote_type,A),rows,A...)
macro nocat(x)
    ex = expand(x)
    if ex.args[1] == :vcat
        ex.args[1] = :(ApproxFun.vnocat)
    elseif ex.args[1] == :hcat
        ex.args[1] = :(ApproxFun.hnocat)
    else
        @assert ex.args[1] == :hvcat
        ex.args[1] = :(ApproxFun.hvnocat)
    end
    esc(ex)
end



## Dynamic functions

struct DFunction <: Function
    f
end
(f::DFunction)(args...) = f.f(args...)

hasnumargs(f::DFunction, k) = hasnumargs(f.f, k)

dynamic(f) = f
dynamic(f::Function) = DFunction(f) # Assume f has to compile every time
