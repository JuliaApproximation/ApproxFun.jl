export clenshaw, @clenshaw

type ClenshawPlan{S,T}
    sp::S
    bk::Vector{T}
    bk1::Vector{T}
    bk2::Vector{T}
    A::Vector{T}
    B::Vector{T}
    C::Vector{T}
end

function ClenshawPlan{S,T}(::Type{T},sp::S,N::Int,n::Int)
    A = T[recA(T,sp,k) for k=0:N-1]
    B = T[recB(T,sp,k) for k=0:N-1]
    C = T[recC(T,sp,k) for k=1:N]
    ClenshawPlan(sp,Array(T,n),Array(T,n),Array(T,n),A,B,C)
end

macro clenshaw(x, c...)
    bk1,bk2 = :(zero(t)),:(zero(t))
    N = length(c)
    for k = N:-1:2
        bk2, bk1 = bk1, :(muladd(t,$bk1,$(esc(c[k]))-$bk2))
    end
    ex = :(muladd(t/2,$bk1,$(esc(c[1]))-$bk2))
    Expr(:block, :(t = $(esc(2))*$(esc(x))), ex)
end

clenshaw{S,T<:Number,U<:Number}(sp::S,c::AbstractVector{T},x::AbstractArray{U}) = clenshaw(c,x,ClenshawPlan(promote_type(T,U),sp,length(c),length(x)))
clenshaw{S,T<:Number,U<:Number}(sp::S,c::AbstractMatrix{T},x::AbstractArray{U}) = clenshaw(c,x,ClenshawPlan(promote_type(T,U),sp,size(c,1),length(x)))
clenshaw{S,T<:Number,U<:Number}(sp::S,c::AbstractMatrix{T},x::U) = clenshaw(c,x,ClenshawPlan(promote_type(T,U),sp,size(c,1),size(c,2)))

clenshaw{S,T<:Number,U<:Number,V<:Number}(c::AbstractVecOrMat{T},x::AbstractArray{U},plan::ClenshawPlan{S,V}) = reshape(clenshaw(c,vec(x),plan),size(x))
clenshaw{S,T<:Number,U<:Number,V<:Number}(c::AbstractVecOrMat{T},x::U,plan::ClenshawPlan{S,V}) = reshape(clenshaw(c,x,plan),size(c,2))

clenshaw!{S,T,U}(sp::S,c::AbstractVector{T},x::AbstractVector{U})=clenshaw!(c,x,ClenshawPlan(promote_type(T,U),sp,length(x)))

function clenshaw{S,U<:Number,V<:Number}(sp::S,c::AbstractVector{U},x::V)
    N,T = length(c),promote_type(U,V)
    if isempty(c)
        return zero(x)
    end

    bk1,bk2 = zero(T),zero(T)
    A,B,C = recA(T,sp,N-1),recB(T,sp,N-1),recC(T,sp,N)
    for k = N:-1:2
        bk2, bk1 = bk1, muladd(muladd(A,x,B),bk1,muladd(-C,bk2,c[k])) # muladd(-C,bk2,muladd(muladd(A,x,B),bk1,c[k])) # (A*x+B)*bk1+c[k]-C*bk2
        A,B,C = recA(T,sp,k-2),recB(T,sp,k-2),recC(T,sp,k-1)
    end
    muladd(muladd(A,x,B),bk1,muladd(-C,bk2,c[1])) # muladd(-C,bk2,muladd(muladd(A,x,B),bk1,c[1])) # (A*x+B)*bk1+c[1]-C*bk2
end

function clenshaw{S,T<:Number,U<:Number,V<:Number}(c::AbstractVector{T},x::AbstractVector{U},plan::ClenshawPlan{S,V})
    N,n = length(c),length(x)
    if isempty(c)
        return zeros(x)
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

function clenshaw{S,T<:Number,U<:Number,V<:Number}(c::AbstractMatrix{T},x::U,plan::ClenshawPlan{S,V})
    N,n = size(c)
    if isempty(c)
        return zeros(U,n)
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

function clenshaw{S,T<:Number,U<:Number,V<:Number}(c::AbstractMatrix{T},x::AbstractVector{U},plan::ClenshawPlan{S,V})
    N,n = size(c)
    @assert n == length(x)
    if isempty(c)
        return zeros(x)
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

function clenshaw!{S,T<:Number,U<:Number,V<:Number}(c::Vector{T},x::Vector{U},plan::ClenshawPlan{S,V})
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


function sineshaw(c::Vector,θ::Number)
    if isempty(c)
        return zero(θ)
    end

    x = 2cos(θ)
    bk1,bk2 = zero(x),zero(x)
    for k = length(c):-1:1
        bk2, bk1 = bk1, muladd(x,bk1,c[k]-bk2)
    end

    sin(θ)*bk1
end
sineshaw(c::Vector,θ::Vector) = promote_type(eltype(c),eltype(θ))[sineshaw(c,θ[k]) for k=1:length(θ)]
sineshaw(c::Vector,θ::Matrix) = promote_type(eltype(c),eltype(θ))[sineshaw(c,θ[k,j]) for k=1:size(θ,1),j=1:size(θ,2)]
