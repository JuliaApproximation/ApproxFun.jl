type ClenshawPlan{S,T}
    sp::S
    bk::Vector{T}
    bk1::Vector{T}
    bk2::Vector{T}
    rα::Vector{T}
    rβ::Vector{T}
    rγ::Vector{T}
end

function ClenshawPlan{S,T}(::Type{T},sp::S,N::Int,n::Int)
    rα = T[recα(T,sp,k) for k=1:N]
    rβ = T[recβ(T,sp,k) for k=1:N+1]
    rγ = T[recγ(T,sp,k) for k=1:N]
    ClenshawPlan(sp,Array(T,n),Array(T,n),Array(T,n),rα,rβ,rγ)
end

function clenshaw{S,U,V}(sp::S,c::AbstractVector{U},x::V)
    N,T = length(c),promote_type(U,V)
    if isempty(c)
        return zero(x)
    end

    bk1,bk2 = zero(T),zero(T)
    rα1,rβ1,rβ2,rγ2 = recα(T,sp,N),recβ(T,sp,N),recβ(T,sp,N+1),recγ(T,sp,N+1)
    for k = N:-1:2
        bk2, bk1 = bk1, muladd((x-rα1)/rβ1,bk1,muladd(-rγ2/rβ2,bk2,c[k]))
        ra1,rβ1,rβ2,rγ2 = recα(T,sp,k-1),recβ(T,sp,k-1),rβ1,recγ(T,sp,k)
    end

    muladd((x-rα1)/rβ1,bk1,muladd(-rγ2/rβ2,bk2,c[1]))
end

clenshaw{S,T,U}(sp::S,c::AbstractVector{T},x::AbstractVector{U}) = clenshaw(c,x,ClenshawPlan(promote_type(T,U),sp,length(c),length(x)))

function clenshaw{S,T,U,V}(c::AbstractVector{T},x::AbstractVector{U},plan::ClenshawPlan{S,V})
    N,n = length(c),length(x)
    if isempty(c)
        return zeros(x)
    end

    bk=plan.bk
    bk1=plan.bk1
    bk2=plan.bk2
    rα=plan.rα
    rβ=plan.rβ
    rγ=plan.rγ

    @inbounds for i = 1:n
        bk1[i] = zero(V)
        bk2[i] = zero(V)
    end

    @inbounds for k = N:-1:2
        ck,rα1,rβ1,rβ2,rγ2 = c[k],rα[k],rβ[k],rβ[k+1],rγ[k+1]
        for i = 1:n
            bk[i] = muladd((x[i]-rα1)/rβ1,bk1[i],muladd(-rγ2/rβ2,bk2[i],ck))
        end
        bk2, bk1, bk = bk1, bk, bk2
    end

    ck,rα1,rβ1,rβ2,rγ2 = c[1],rα[1],rβ[1],rβ[2],rγ[2]
    @inbounds for i = 1:n
        bk[i] = muladd((x[i]-rα1)/rβ1,bk1[i],muladd(-rγ2/rβ2,bk2[i],ck))
    end

    bk
end


########################################################################################################


#Clenshaw routine for many Funs, x is a vector of same number of funs
#each fun is a column
clenshaw{T<:Number}(c::Array{T,2},x::Vector{T})=clenshaw(c,x,ClenshawPlan(T,size(c)[2]))
function clenshaw{T<:Number}(c::Array{T,2},x::Vector{T},plan::ClenshawPlan{T})
    bk=plan.bk
    bk1=plan.bk1
    bk2=plan.bk2


    m,n=size(c) # m is # of coefficients, n is # of funs

    for i = 1:n
        @inbounds bk1[i] = zero(T)
        @inbounds bk2[i] = zero(T)
        @inbounds bk[i] = zero(T)
    end

    for k=m:-1:2
        for j=1:n
            ck = c[k,j]

            @inbounds bk[j] = ck + 2x[j]*bk1[j] - bk2[j]
        end

        bk2, bk1, bk = bk1, bk, bk2
    end


    for j = 1:n
        ce = c[1,j]
        @inbounds bk[j] = ce + x[j]*bk1[j] - bk2[j]
    end

    bk
end


#Clenshaw routine for many Funs
#each fun is a column
clenshaw{T<:Number}(c::Array{T,2},x::Number)=clenshaw(c,x,ClenshawPlan(promote_type(T,typeof(x)),size(c)[2]))
function clenshaw{T<:Number}(c::Array{T,2},x::Number,plan::ClenshawPlan{T})
    bk=plan.bk
    bk1=plan.bk1
    bk2=plan.bk2

    n=size(c)[2] #number of funs
    m=size(c)[1] #number of coefficients

    for i = 1:n
        @inbounds bk1[i] = zero(T)
        @inbounds bk2[i] = zero(T)
        @inbounds bk[i] = zero(T)
    end

    for k=m:-1:2
        for j=1:n
            @inbounds ck = c[k,j]

            @inbounds bk[j] = ck + 2x * bk1[j] - bk2[j]
        end

        bk2, bk1, bk = bk1, bk, bk2
    end

    for j = 1:n
        ce = c[1,j]
        @inbounds bk[j] = ce + x * bk1[j] - bk2[j]
    end

    bk
end



# overwrite x
clenshaw!{T<:Number,M<:Number}(c::Vector{T},x::Vector{M})=clenshaw!(c,x,ClenshawPlan(promote_type(T,M),length(x)))
function clenshaw!{T<:Number,M<:Number,Q<:Number}(c::Vector{T},x::Vector{M},plan::ClenshawPlan{Q})
    n = length(x)

    if isempty(c)
        for k=1:n
            x[k]=zero(T)
        end
        return x
    end

    bk=plan.bk
    bk1=plan.bk1
    bk2=plan.bk2

#    x=2x
    for i = 1:n
        @inbounds bk1[i] = zero(Q)
        @inbounds bk2[i] = zero(Q)
        @inbounds bk[i] = zero(Q)
    end

    for k in  length(c):-1:2
        ck = c[k]
        for i in 1 : n
            @inbounds bk[i] = ck + 2x[i] * bk1[i] - bk2[i]
        end
        bk2, bk1, bk = bk1, bk, bk2
    end

    ce = c[1]
    for i in 1 : n
        @inbounds  x[i] = ce + x[i] * bk1[i] - bk2[i]
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
