
import Base.chop

# Used for spaces not defined yet
immutable UnsetNumber <: Number  end
Base.promote_rule{N<:Number}(::Type{UnsetNumber},::Type{N})=N

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
eps{T<:Number}(z::Type{Vector{T}})=eps(T)

dotu(f::Vector{Complex{Float64}},g::Vector{Complex{Float64}})=BLAS.dotu(f,g)
dotu{N<:Real}(f::Vector{Complex{Float64}},g::Vector{N})=dot(conj(f),g)
dotu{N<:Real,T<:Number}(f::Vector{N},g::Vector{T})=dot(f,g)


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


function interlace{T<:Number}(v::Vector{Vector{T}})
    n=length(v)
    l=mapreduce(length,max,v)
    ret=zeros(T,n*l)

    for k=1:n
        ret[k:n:k+(length(v[k])-1)*n]=v[k]
    end
    ret
end

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

function interlace(a::Vector,b::Vector)
    na=length(a);nb=length(b)
    T=promote_type(eltype(a),eltype(b))
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




## svfft

##FFT That interlaces coefficients

plan_svfft(x::Vector) = wrap_fft_plan(plan_fft(x))
plan_isvfft(x::Vector) = wrap_fft_plan(plan_ifft(x))

function svfft(v::Vector,plan)
    n=length(v)
    v=plan(v)/n
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

    plan(n*v)
end





## tridiagonal ql

function tridql!(L::Matrix)
    n=size(L,1)

  # Now we do QL for the compact part in the top left
    Q = eye(eltype(L),n)
    for i = n:-1:2
        nrm=sqrt(L[i-1,i]^2+L[i,i]^2)
        c,s = L[i,i]/nrm, L[i-1,i]/nrm
        if i > 2
            L[i-1:i,i-2:i] = [c -s; s c]*L[i-1:i,i-2:i]
            L[i-1,i]=0
        else
            L[i-1:i,i-1:i] = [c -s; s c]*L[i-1:i,i-1:i]
            L[i-1,i]=0
        end
        G=eye(eltype(L),n)
        G[i-1:i,i-1:i]=[c s; -s c]
        Q = Q*G
    end
    Q,L
end
