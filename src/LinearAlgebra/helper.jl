export dotu
import Base.chop

# Used for spaces not defined yet
immutable UnsetNumber <: Number  end
Base.promote_rule{N<:Number}(::Type{UnsetNumber},::Type{N}) = N

# Test the number of arguments a function takes
hasnumargs(f,k) = applicable(f,zeros(k)...)


isapprox(a...;kwds...) = Base.isapprox(a...;kwds...)
isapprox(a::Vec,b::Vec;kwds...) = isapprox([a...],[b...];kwds...)

# fast implementation of isapprox with atol a non-keyword argument in most cases
isapprox_atol(a,b,atol;kwds...) = isapprox(a,b;atol=atol,kwds...)
isapprox_atol(a::Vec,b::Vec,atol::Real=0;kwds...) = isapprox_atol([a...],[b...],atol;kwds...)
function isapprox_atol(x::Number, y::Number, atol::Real=0; rtol::Real=Base.rtoldefault(x,y))
    x == y || (isfinite(x) && isfinite(y) && abs(x-y) <= atol + rtol*max(abs(x), abs(y)))
end
function isapprox_atol{T<:Number,S<:Number}(x::AbstractArray{T}, y::AbstractArray{S},atol::Real=0; rtol::Real=Base.rtoldefault(T,S), norm::Function=vecnorm)
    d = norm(x - y)
    if isfinite(d)
        return d <= atol + rtol*max(norm(x), norm(y))
    else
        # Fall back to a component-wise approximate comparison
        return all(ab -> isapprox(ab[1], ab[2]; rtol=rtol, atol=atol), zip(x, y))
    end
end


# This creates ApproxFun.real, ApproxFun.eps and ApproxFun.dou
# which we override for default julia types
real(x...) = Base.real(x...)
real(::Type{UnsetNumber}) = UnsetNumber
real{T<:Real}(::Type{T}) = T
real{T<:Real}(::Type{Complex{T}}) = T
real{T<:Real,n}(::Type{Array{T,n}}) = Array{T,n}
real{T<:Complex,n}(::Type{Array{T,n}}) = Array{real(T),n}
real{N,T<:Real}(::Type{Vec{N,T}}) = Vec{N,T}
real{N,T<:Complex}(::Type{Vec{N,T}}) = Vec{N,real(T)}


eps(x...) = Base.eps(x...)
eps{T<:Real}(::Type{Complex{T}}) = eps(real(T))
eps{T<:Real}(z::Complex{T}) = eps(abs(z))
eps{T<:Real}(::Type{Dual{Complex{T}}}) = eps(real(T))
eps{T<:Real}(z::Dual{Complex{T}}) = eps(abs(z))
eps{T<:Number}(::Type{Vector{T}}) = eps(T)
eps{k,T<:Number}(::Type{Vec{k,T}}) = eps(T)


isnan(x) = Base.isnan(x)
isnan(x::Vec) = map(isnan,x)


# BLAS


# implement muladd default
muladd(a,b,c) = a*b+c
muladd(a::Number,b::Number,c::Number) = Base.muladd(a,b,c)


for TYP in (:Float64,:Float32,:Complex128,:Complex64)
    @eval scal!{T<:$TYP}(n::Integer,cst::$TYP,ret::DenseArray{T},k::Integer) =
            BLAS.scal!(n,cst,ret,k)
end


scal!{T<:BlasFloat}(n::Integer,cst::BlasFloat,ret::DenseArray{T},k::Integer) =
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

function reverseeven!(x::Vector)
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

function negateeven!(x::Vector)
    @inbounds @simd for k = 2:2:length(x)
        x[k] *= -1
    end
    x
end

#checkerboard, same as applying negativeeven! to all rows then all columns
function negateeven!(X::Matrix)
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

alternatesign(v::Vector)=alternatesign!(copy(v))

alternatingvector(n::Integer) = 2*mod([1:n],2) .- 1

function alternatingsum(v::Vector)
    ret = zero(eltype(v))
    s = 1
    @inbounds for k=1:length(v)
        ret+=s*v[k]
        s*=-1
    end

    ret
end

# Sum Hadamard product of vectors up to minimum over lengths
function mindotu(a::Vector,b::Vector)
    ret,m = zero(promote_type(eltype(a),eltype(b))),min(length(a),length(b))
    @inbounds @simd for i=m:-1:1 ret += a[i]*b[i] end
    ret
end


# efficiently resize a Matrix.  Note it doesn't change the input ptr
function unsafe_resize!(W::Matrix,::Colon,m::Integer)
    if m == size(W,2)
        W
    else
        n=size(W,1)
        reshape(resize!(vec(W),n*m),n,m)
    end
end

function unsafe_resize!(W::Matrix,n::Integer,::Colon)
    N=size(W,1)
    if n == N
        W
    elseif n < N
        W[1:n,:]
    else
        m=size(W,2)
        ret=Array(eltype(W),n,m)
        ret[1:N,:] = W
        ret
    end
end

function unsafe_resize!(W::Matrix,n::Integer,m::Integer)
    N=size(W,1)
    if n == N
        unsafe_resize!(W,:,m)
    else
        unsafe_resize!(unsafe_resize!(W,n,:),:,m)
    end
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
	   ret=Vector{T}(n)
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
	   Matrix{T}(n,m)  #fixes weird julia bug when T==None
    else
        ret = Matrix{T}(n,m)
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
    b=Vector{Vector{T}}(length(v))
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

function cycle_rotate!(v::Vector, leader::Int, it::Int, twom::Int)
    i = nextindex(leader, twom)
    while i != leader
        idx1, idx2 = it + i - 1, it + leader - 1
        @inbounds v[idx1], v[idx2] = v[idx2], v[idx1]
        i = nextindex(i, twom)
    end
    v
end

function right_cyclic_shift!(v::Vector, it::Int, m::Int, n::Int)
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
function interlace!(v::Vector,offset::Int)
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

function slnorm(u::AbstractMatrix,r::Range,::Colon)
    ret = zero(real(eltype(u)))
    for k=r
        @simd for j=1:size(u,2)
            #@inbounds
            ret=max(norm(u[k,j]),ret)
        end
    end
    ret
end


function slnorm(m::AbstractMatrix,kr::Range,jr::Range)
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

slnorm(m::Matrix,kr::Range,jr::Integer) = slnorm(m,kr,jr:jr)
slnorm(m::Matrix,kr::Integer,jr::Range) = slnorm(m,kr:kr,jr)


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


slnorm(m::AbstractMatrix,k::Integer,::Colon) = slnorm(m,k,1:size(m,2))
slnorm(m::AbstractMatrix,::Colon,j::Integer) = slnorm(m,1:size(m,1),j)



## New Inf

# angle is π*a where a is (false==0) and (true==1)
immutable Infinity{T} <: Number
    angle::T
end

Infinity() = Infinity(false)
const ∞ = Infinity()


Base.isinf(::Infinity) = true
Base.isfinite(::Infinity) = false
Base.sign{B<:Integer}(y::Infinity{B}) = mod(y.angle,2)==0?1:-1
Base.angle(x::Infinity) = π*x.angle

Base.zero{B}(::Infinity{B}) = zero(B)
Base.one{B}(::Infinity{B}) = one(B)


function Base.show{B<:Integer}(io::IO, y::Infinity{B})
    if sign(y) == 1
        print(io, "∞")
    else
        print(io, "-∞")
    end
end

Base.show(io::IO,x::Infinity) = print(io,"$(exp(im*π*x.angle))∞")

Base.promote_rule{R<:Number,B}(::Type{Infinity{B}},::Type{R}) = Number

==(x::Infinity,y::Infinity) = x.angle == y.angle
for TYP in (:Dual,:Number)
    @eval begin
        ==(x::Infinity,y::$TYP) = isinf(y) && angle(y) == angle(x)
        ==(y::$TYP,x::Infinity) = x == y
    end
end
Base.isless(x::Infinity{Bool},y::Infinity{Bool}) = x.angle && !y.angle
Base.isless(x::Number,y::Infinity{Bool}) = y.angle && (x ≠ ∞)
Base.isless(x::Infinity{Bool},y::Number) = !x.angle && (y ≠ -∞)


-{B<:Integer}(y::Infinity{B}) = sign(y)==1?Infinity(one(B)):Infinity(zero(B))

function +{B}(x::Infinity{B},y::Infinity{B})
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


# $ is xor
*(a::Infinity{Bool},b::Infinity{Bool}) = Infinity(a.angle $ b.angle)
*(a::Infinity,b::Infinity) = Infinity(a.angle + b.angle)

for T in (:Dual,:Bool,:Integer,:AbstractFloat)
    @eval begin
        *(a::$T,y::Infinity) = a>0?y:(-y)
        *(y::Infinity,a::$T) = a*y
    end
end

*(a::Number,y::Infinity) = Infinity(y.angle+angle(a)/π)
*(y::Infinity,a::Number) = a*y

for OP in (:fld,:cld,:div)
  @eval Base.$OP(y::Infinity,a::Number) = y*(1/sign(a))
end

Base.min{B<:Integer}(x::Infinity{B},y::Infinity{B}) = sign(x)==-1?x:y
Base.max{B<:Integer}(x::Infinity{B},::Infinity{B}) = sign(x)==1?x:y
Base.min{B<:Integer}(x::Real,y::Infinity{B}) = sign(y)==1?x:y
Base.min{B<:Integer}(x::Infinity{B},y::Real) = min(y,x)
Base.max{B<:Integer}(x::Real,y::Infinity{B}) = sign(y)==1?y:x
Base.max{B<:Integer}(x::Infinity{B},y::Real) = max(y,x)

for OP in (:<,:<=)
    @eval begin
        $OP{B<:Integer}(x::Real,y::Infinity{B}) = sign(y)==1
        $OP{B<:Integer}(y::Infinity{B},x::Real) = sign(y)==-1
    end
end

for OP in (:>,:>=)
    @eval begin
        $OP{B<:Integer}(x::Real,y::Infinity{B}) = sign(y)==-1
        $OP{B<:Integer}(y::Infinity{B},x::Real) = sign(y)==1
    end
end


# Take -- iterate through the first n elements

immutable Take{I,T} <: AbstractVector{T}
    xs::I
    n::Int
    function Take(xs,n)
        @assert n ≤ length(xs)
        new(xs,n)
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

eltype{I}(::Type{Take{I}}) = eltype(I)
Base.length(t::Take) = t.n
Base.size(t::Take) = (length(t),)
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

function Base.sum(it::Take)
    ret = zero(eltype(it))
    for a in it
        ret += a
    end
    ret
end

pad(it::Take,n::Integer) = pad!(collect(it),n)



# Re-implementation of Base iterators
# to use ∞ and allow getindex

abstract AbstractRepeated{T}

Base.eltype{T}(::Type{AbstractRepeated{T}}) = T
Base.eltype{R<:AbstractRepeated}(::Type{R}) = eltype(super(R))
Base.eltype{T}(::AbstractRepeated{T}) = T

Base.step(::AbstractRepeated) = 0

Base.start(::AbstractRepeated) = nothing
Base.next(it::AbstractRepeated,state) = value(it),nothing
Base.done(::AbstractRepeated,state) = false

Base.length(::AbstractRepeated) = ∞

getindex(it::AbstractRepeated,k::Integer) = value(it)
getindex(it::AbstractRepeated,k::Range) = take(it,length(k))


Base.maximum(r::AbstractRepeated) = value(r)
Base.minimum(r::AbstractRepeated) = value(r)

Base.sum{AR<:AbstractRepeated}(it::Take{AR}) = it.n*value(it.xs)

immutable ZeroRepeated{T} <: AbstractRepeated{T} end

ZeroRepeated{T}(::Type{T}) = ZeroRepeated{T}()

value{T}(::ZeroRepeated{T}) = zero(T)
Base.sum(r::ZeroRepeated) = value(r)

immutable Repeated{T} <: AbstractRepeated{T}
    x::T
    function Repeated(x::T)
        # TODO: Add ZeroRepeated type.
        if x == zero(T)
            error("Zero repeated not supported to maintain type stability")
        end

        new(x)
    end
end

Repeated(x) = Repeated{typeof(x)}(x)


value(r::Repeated) = r.x

Base.sum(r::Repeated) = r.x > 0 ? ∞ : -∞



function repeated(x)
    if x == zero(x)
        error("Use ZeroRepeated to repeat zeros")
    end
    Repeated(x)
end
repeated(x,::Infinity{Bool}) = repeated(x)
repeated(x,m::Integer) = take(repeated(x),m)


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
Base.eltype{AS<:AbstractCount}(::Type{AS}) = eltype(supertype(AS))
Base.eltype{S}(::AbstractCount{S}) = S

Base.step(it::Count) = it.step
Base.step(it::UnitCount) = 1

Base.start(it::AbstractCount) = it.start
Base.next(it::AbstractCount, state) = (state, state + step(it))
Base.done(it::AbstractCount, state) = false

Base.length(it::AbstractCount) = ∞

getindex(it::Count,k) = it.start + it.step*(k-1)
getindex(it::UnitCount,k) = (it.start-1) + k
getindex(it::AbstractRepeated,k::AbstractCount) = it

# use reindex, copied from Base
reindex(V, idxs::Tuple{AbstractCount, Vararg{Any}}, subidxs::Tuple{Any, Vararg{Any}}) =
    (Base.@_propagate_inbounds_meta; (idxs[1][subidxs[1]], reindex(V, tail(idxs), tail(subidxs))...))

# special functions
Base.maximum{S<:Real}(it::UnitCount{S}) = ∞
Base.minimum{S<:Real}(it::UnitCount{S}) = it.start

Base.maximum{S<:Real}(it::Count{S}) = it.step > 0 ? ∞ : it.start
Base.minimum{S<:Real}(it::Count{S}) = it.step < 0 ? -∞ : it.start

Base.sum(it::UnitCount) = ∞
Base.sum(it::Count) = it.step > 0 ? ∞ : -∞

Base.last{S<:Real}(it::UnitCount{S}) = ∞
Base.last{S<:Real}(it::Count{S}) =
    it.step > 0 ? ∞ : (it.step < 0 ? -∞ : error("zero step not supported"))


function Base.colon(a::Real,b::Infinity{Bool})
    if b.angle
        throw(ArgumentError("Cannot create $a:-∞"))
    end
    countfrom(a)
end
function Base.colon(a::Infinity{Bool},st::AbstractFloat,b::Infinity{Bool})
    if a ≠ b
        throw(ArgumentError("Cannot create $a:$st:$b"))
    end
    [a]
end
function Base.colon(a::Real,st::Real,b::Infinity{Bool})
    if st == 0
        throw(ArgumentError("step cannot be zero"))
    elseif b.angle == st > 0
        throw(ArgumentError("Cannot create $a:$st:$b"))
    else
        countfrom(a,st)
    end
end




immutable CumSumIterator{CC}
    iterator::CC
end

Base.eltype{S}(::Type{CumSumIterator{S}}) = eltype(S)
Base.eltype(CC::CumSumIterator) = eltype(CC.iterator)

Base.start(it::CumSumIterator) = (0,start(it.iterator))
function Base.next(it::CumSumIterator, state)
    a,nx_st=next(it.iterator,state[2])
    cs=state[1]+a
    (cs,(cs,nx_st))
end
Base.done(it::CumSumIterator, state) = done(it.iterator,state[2])

Base.length(it::CumSumIterator) = length(it.iterator)

getindex{AC<:UnitCount}(it::CumSumIterator{AC},k) = it.iterator.start*k + ((k*(k-1))÷2)
getindex{AC<:Count}(it::CumSumIterator{AC},k) = it.iterator.start*k + step(it.iterator)*((k*(k-1))÷2)





## BandedMatrix



pad!(A::BandedMatrix,n,::Colon) = pad!(A,n,n+A.u)  # Default is to get all columns
columnrange(A,row::Integer) = max(1,row+bandinds(A,1)):row+bandinds(A,2)



## Store iterator
type CachedIterator{T,IT,ST}
    iterator::IT
    storage::Vector{T}
    state::ST
    length::Int

    CachedIterator(it::IT) = new(it,Vector{T}(),start(it),0)
end

CachedIterator(it) = CachedIterator{eltype(it),typeof(it),typeof(start(it))}(it)

function Base.resize!(it::CachedIterator,n::Integer)
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


Base.eltype{T}(it::CachedIterator{T}) = T
Base.start(it::CachedIterator) = 1
Base.next(it::CachedIterator,st::Int) = (it[st],st+1)
Base.done(it::CachedIterator,st::Int) = st == it.length + 1 &&
                                        done(it.iterator,it.state)

function getindex(it::CachedIterator,k)
    mx = maximum(k)
    if mx > length(it) || mx < 1
        throw(BoundsError(it,k))
    end
    resize!(it,isempty(k)?0:mx).storage[k]
end
function Base.findfirst(f::Function,A::CachedIterator)
    k=1
    for c in A
        if f(c)
            return k
        end
        k+=1
    end
    return 0
end

function Base.findfirst(A::CachedIterator,x)
    k=1
    for c in A
        if c == x
            return k
        end
        k+=1
    end
    return 0
end

Base.length(A::CachedIterator) = length(A.iterator)


# The following don't need caching
cache{T<:Number}(A::Vector{T}) = A
cache(A::Range) = A
cache(A::AbstractCount) = A



## From Julia v0.5 code


# flatten an iterator of iterators
# we add indexing

immutable Flatten{I}
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

Base.eltype(f::Flatten) = mapreduce(eltype,promote_type,f.it)


Base.start(f::Flatten) = 1, map(start,f.it)

function Base.next(f::Flatten, state)
    k, sts = state
    if !done(f.it[k],sts[k])
        a,nst = next(f.it[k],sts[k])
        a, (k,(sts[1:k-1]...,nst,sts[k+1:end]...))
    else
        next(f,(k+1,sts))
    end
end


@inline Base.done(f::Flatten, state) =
    state[1] == length(f) && done(f.it[end],state[2][end])

Base.length(f::Flatten) = mapreduce(length,+,f.it)

Base.eachindex(f::Flatten) = 1:length(f)

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

Base.sum(f::Flatten) = mapreduce(sum,+,f.it)


Base.maximum(f::Flatten) = mapreduce(maximum,max,f.it)
Base.minimum(f::Flatten) = mapreduce(minimum,min,f.it)


## Iterator Algebra

broadcast(op,f::Flatten,c...) = Flatten(map(it->op(it,c...),f.it))

if VERSION < v"0.6.0-dev"
    for TYP in (:Flatten, :AbstractRepeated), op in (:+,:-,:*,:/)
        dop = parse("."*string(op))
        @eval begin
            $dop(a::$TYP,b::$TYP) = broadcst($op,a,b)
            $dop(f::$TYP,c::Number) = broadcast($op,f,c)
            $dop(c::Number,f::$TYP) = broadcast($op,c,f)
        end
    end
end


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

typealias InfiniteIterators Union{AbstractRepeated,AbstractCount,Flatten}

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


Base.cumsum(r::Repeated) = r.x:r.x:(r.x>0?∞:-∞)
Base.cumsum(r::Repeated{Bool}) = 1:∞
Base.cumsum(r::ZeroRepeated) = r
Base.cumsum(r::AbstractCount) = CumSumIterator(r)




function Base.cumsum(f::Flatten)
    cs=zero(eltype(f))
    its = Vector{eltype(f.it)}(0)
    for it in f.it[1:end-1]
        c=cumsum(cs+it)
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
