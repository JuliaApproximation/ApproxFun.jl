

##
# Represent an m x n rank-r matrix
# A = U V^T
##


balance!(U::Matrix{T},V::Matrix{T},m::Int,n::Int,r::Int) where {T<:Union{Integer,Rational}} = U,V
function balance!(U::Matrix{T},V::Matrix{T},m::Int,n::Int,r::Int) where T
    for k=1:r
        uk = zero(T)
        for i=1:m
            @inbounds uk += abs2(U[i,k])
        end
        vk = zero(T)
        for j=1:n
            @inbounds vk += abs2(V[j,k])
        end
        uk,vk = sqrt(uk),sqrt(vk)
        σk = sqrt(uk*vk)
        if abs2(uk) ≥ eps(T)^2 && abs2(vk) ≥ eps(T)^2
            uk,vk = σk/uk,σk/vk
            for i=1:m
                @inbounds U[i,k] *= uk
            end
            for j=1:n
                @inbounds V[j,k] *= vk
            end
        end
    end
    U,V
end

function refactorsvd!(U::Matrix{S},Σ::Vector{T},V::Matrix{S}) where {S,T}
    conj!(V)
    σmax = Σ[1]
    r=max(1,count(s->s>10σmax*eps(T),Σ))
    m,n = size(U,1),size(V,1)
    for k=1:r
        σk = sqrt(Σ[k])
        for i=1:m
            @inbounds U[i,k] *= σk
        end
        for j=1:n
            @inbounds V[j,k] *= σk
        end
    end
    r
end

# constructors

function pad!(L::LowRankMatrix,n::Integer,::Colon)
    L.U=pad(L.U,n,:)
    L
end
function pad!(L::LowRankMatrix,::Colon,m::Integer)
    L.V=pad(L.V,m,:)
    L
end
pad!(L::LowRankMatrix,n::Integer,m::Integer) = pad!(pad!(L,n,:),:,m)
