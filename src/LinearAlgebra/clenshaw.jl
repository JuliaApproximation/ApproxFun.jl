export clenshaw, @clenshaw


##
# ClenshawPlan can be used for multiple evaluations with different functions
# that have the same length
##

mutable struct ClenshawPlan{S,T}
    sp::S
    bk::Vector{T}
    bk1::Vector{T}
    bk2::Vector{T}
    A::Vector{T}
    B::Vector{T}
    C::Vector{T}
end

function ClenshawPlan(::Type{T},sp,N::Int,n::Int) where T
    A = T[recA(T,sp,k) for k=0:N-1]
    B = T[recB(T,sp,k) for k=0:N-1]
    C = T[recC(T,sp,k) for k=1:N]
    ClenshawPlan(sp,Array{T}(undef,n),Array{T}(undef,n),Array{T}(undef,n),A,B,C)
end

macro clenshaw(x, c...)
    a, b = :(zero(t)), :(zero(t))
    as = []
    N = length(c)
    for k = N:-1:2
        ak = Symbol("a",k)
        push!(as, :($ak = $a))
        a = :(muladd(t,$a,$(esc(c[k]))-$b))
        b = :($ak)
    end
    ex = Expr(:block,as...,:(muladd(t/2,$a,$(esc(c[1]))-$b)))
    Expr(:block, :(t = $(esc(2))*$(esc(x))), ex)
end

clenshaw(x,c) = clenshaw_halved(2*x, c)

@generated function clenshaw_halved(x, c::StaticArrays.SVector{N,R}) where {N, R}
    a, b = :(zero(R)), :(zero(R))
    as = []
    for k = N:-1:2
        ak = Symbol("a",k)
        push!(as, :($ak = $a))
        a = :(muladd(x,$a,c[$k]-$b))
        b = :($ak)
    end
    Expr(:block,
    as...,
    :(muladd(x/2,$a,c[1]-$b)))
end

for TYP in (:AbstractVector,:AbstractMatrix)
    @eval clenshaw(c::$TYP,x::AbstractArray,plan::ClenshawPlan) = reshape(clenshaw(c,vec(x),plan),size(x))
end


function clenshaw(c::AbstractVector,x::AbstractVector,plan::ClenshawPlan{S,V}) where {S,V}
    N,n = length(c),length(x)
    if isempty(c)
        return zeros(V,length(x))
    end

    bk=plan.bk
    bk1=plan.bk1
    bk2=plan.bk2
    A=plan.A
    B=plan.B
    C=plan.C

    @inbounds for i = 1:n
        bk1[i] = zero(V)
        bk2[i] = zero(V)
    end

    @inbounds for k = N:-1:2
        ck,Ak,Bk,Ck = c[k],A[k],B[k],C[k]
        for i = 1:n
            bk[i] = muladd(muladd(Ak,x[i],Bk),bk1[i],muladd(-Ck,bk2[i],ck))
        end
        bk2, bk1, bk = bk1, bk, bk2
    end

    ck,Ak,Bk,Ck = c[1],A[1],B[1],C[1]
    @inbounds for i = 1:n
        bk[i] = muladd(muladd(Ak,x[i],Bk),bk1[i],muladd(-Ck,bk2[i],ck))
    end

    bk
end


#Clenshaw routine for many Funs, x is a vector of same number of funs
#each fun is a column

function clenshaw(c::AbstractMatrix,x,plan::ClenshawPlan{S,V}) where {S,V}
    N,n = size(c)
    if isempty(c)
        return zeros(V,n)
    end

    bk=plan.bk
    bk1=plan.bk1
    bk2=plan.bk2
    A=plan.A
    B=plan.B
    C=plan.C

    @inbounds for i = 1:n
        bk1[i] = zero(V)
        bk2[i] = zero(V)
    end

    @inbounds for k = N:-1:2
        Ak,Bk,Ck = A[k],B[k],C[k]
        for i = 1:n
            cki = c[k,i]
            bk[i] = muladd(muladd(Ak,x,Bk),bk1[i],muladd(-Ck,bk2[i],cki))
        end
        bk2, bk1, bk = bk1, bk, bk2
    end

    Ak,Bk,Ck = A[1],B[1],C[1]
    @inbounds for i = 1:n
        ci = c[1,i]
        bk[i] = muladd(muladd(Ak,x,Bk),bk1[i],muladd(-Ck,bk2[i],ci))
    end

    bk
end

function clenshaw(c::AbstractMatrix,x::AbstractVector,plan::ClenshawPlan{S,V}) where {S,V}
    N,n = size(c)
    @assert n == length(x)
    if isempty(c)
        return zeros(V,length(x))
    end

    bk=plan.bk
    bk1=plan.bk1
    bk2=plan.bk2
    A=plan.A
    B=plan.B
    C=plan.C

    @inbounds for i = 1:n
        bk1[i] = zero(V)
        bk2[i] = zero(V)
    end

    @inbounds for k = N:-1:2
        Ak,Bk,Ck = A[k],B[k],C[k]
        for i = 1:n
            cki = c[k,i]
            bk[i] = muladd(muladd(Ak,x[i],Bk),bk1[i],muladd(-Ck,bk2[i],cki))
        end
        bk2, bk1, bk = bk1, bk, bk2
    end

    Ak,Bk,Ck = A[1],B[1],C[1]
    @inbounds for i = 1:n
        ci = c[1,i]
        bk[i] = muladd(muladd(Ak,x[i],Bk),bk1[i],muladd(-Ck,bk2[i],ci))
    end

    bk
end


# overwrite x

function clenshaw!(c::AbstractVector,x::AbstractVector,plan::ClenshawPlan{S,V}) where {S,V}
    N,n = length(c),length(x)
    if isempty(c)
        for k=1:n
            x[k]=zero(V)
        end
    end

    bk=plan.bk
    bk1=plan.bk1
    bk2=plan.bk2
    A=plan.A
    B=plan.B
    C=plan.C

    @inbounds for i = 1:n
        bk1[i] = zero(V)
        bk2[i] = zero(V)
    end

    @inbounds for k = N:-1:2
        ck,Ak,Bk,Ck = c[k],A[k],B[k],C[k]
        for i = 1:n
            bk[i] = muladd(muladd(Ak,x[i],Bk),bk1[i],muladd(-Ck,bk2[i],ck))
        end
        bk2, bk1, bk = bk1, bk, bk2
    end

    ck,Ak,Bk,Ck = c[1],A[1],B[1],C[1]
    @inbounds for i = 1:n
        x[i] = muladd(muladd(Ak,x[i],Bk),bk1[i],muladd(-Ck,bk2[i],ck))
    end

    x
end


function sineshaw(c::AbstractVector,θ)
    if isempty(c)
        return zero(θ)
    end

    x = 2cos(θ)
    bk1,bk2 = zero(x),zero(x)
    @inbounds for k = length(c):-1:1
        bk2, bk1 = bk1, muladd(x,bk1,c[k]-bk2)
    end

    sin(θ)*bk1
end
sineshaw(c::AbstractVector,θ::AbstractVector) = promote_type(eltype(c),eltype(θ))[sineshaw(c,θ[k]) for k=1:length(θ)]
sineshaw(c::AbstractVector,θ::AbstractMatrix) = promote_type(eltype(c),eltype(θ))[sineshaw(c,θ[k,j]) for k=1:size(θ,1),j=1:size(θ,2)]
