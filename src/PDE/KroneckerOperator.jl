export KroneckerOperator


# gives zero block banded matrix, where the blocks are increasing size
# the bandwidths are assumed to be constant

function blockbandzeros(T,n,m::Integer,Al,Au,Bl,Bu)
    l=Al+Bl;u=Au+Bu    
    ret=BandedMatrix(BandedMatrix{T},n,m,l,u)
    
    for k=1:n,j=max(1,k-l):min(m,k+u)
        nl=min(Al,Bu+k-j);nu=min(Au,Bl+j-k)
        ret[k,j]=bazeros(T,k,j,nl,nu)
    end 
    
    ret
end

blockbandzeros(T,n,m::Colon,Al,Au,Bl,Bu)=blockbandzeros(T,n,n+Au+Bu,Al,Au,Bl,Bu)
blockbandzeros(T,n,m,Alu,Blu)=blockbandzeros(T,n,m,-Alu[1],Alu[2],-Blu[1],Blu[2])


##########
# Convert a block banded matrix to a full matrix
# TODO: Don't assume block banded matrix has i x j blocks
###########

function Base.convert{T,V<:Number}(::Type{Matrix{T}},K::BandedMatrix{BandedMatrix{V}})
    n=size(K,1)
    m=size(K,2)

    ret=zeros(T,div(n*(n+1),2),div(m*(m+1),2))

    for k=1:n,j=max(1,k-K.l):min(m,k+K.u)
        for κ=1:k,ξ=max(1,κ-K[k,j].l):min(j,κ+K[k,j].u)
            ret[div((k-1)*k,2)+κ,div((j-1)*j,2)+ξ]=K[k,j][κ,ξ]
        end
    end 
    ret
end


##############
# BivariateOperator represents a block banded operator 
# the (i,j) block is a i x j BandedMatrix
##############


abstract BivariateOperator{T} <: BandedOperator{BandedMatrix{T}}


##########
# KroneckerOperator gives the kronecker product of two 1D operators
#########

immutable KroneckerOperator{S,V,DS,RS,T}<: BivariateOperator{T}
    ops::(S,V) 
    domainspace::DS
    rangespace::RS
end

KroneckerOperator(A,B,ds,rs)=KroneckerOperator{typeof(A),typeof(B),typeof(ds),typeof(rs),promote_type(eltype(A),eltype(B))}((A,B),ds,rs)
KroneckerOperator(A,B)=KroneckerOperator(A,B,domainspace(A)⊗domainspace(B),rangespace(A)⊗rangespace(B))



bandinds(K::KroneckerOperator)=bandinds(K.ops[1],1)+bandinds(K.ops[2],1),bandinds(K.ops[1],2)+bandinds(K.ops[2],2)
domainspace(K::KroneckerOperator)=K.domainspace
rangespace(K::KroneckerOperator)=K.rangespace

function kronaddentries!(A,B,M,kr::Range)
    m=max(size(A,2),size(B,2))
    l=A.l+B.l;u=A.u+B.u

    for k=kr,j=max(1,k-l):k+u
        nl=min(A.l,B.u+k-j);nu=min(A.u,B.l+j-k)
        for κ=1:k,ξ=max(1,κ-nl):min(j,κ+nu)
            M[k,j][κ,ξ]+=A[κ,ξ]*B[k-κ+1,j-ξ+1]
        end
    end
    M
end

addentries!(K::KroneckerOperator,A,kr::Range)=kronaddentries!(slice(K.ops[1],1:last(kr),:),slice(K.ops[2],1:last(kr),:),A,kr)


BandedMatrix(K::KroneckerOperator,n::Integer)=addentries!(K,blockbandzeros(Float64,n,:,bandinds(K.ops[1]),bandinds(K.ops[2])),1:n)
function BandedMatrix(K::KroneckerOperator,kr::UnitRange,::Colon)
    @assert first(kr)==1
    BandedMatrix(K,last(kr))
end





##########
# Multiply a block banded matrix by a vector, where the vector is assumed to 
# decompose into blocks
# TODO: Don't assume block banded matrix has i x j blocks
###########

function *{T,V<:Number}(M::BandedMatrix{BandedMatrix{T}},v::Vector{V})
    n,m=size(M)
    r=zeros(promote_type(T,V),div(n*(n+1),2))

    for j=1:m-1
        vj=v[fromtensorblock(j)]

        for k=max(1,j-M.u):j+M.l
            r[fromtensorblock(k)]+=M[k,j]*vj
        end
    end

    # pad so block size matches
    vj=pad!(v[first(fromtensorblock(m)):end],m)
    for k=max(1,m-M.u):m+M.l
        r[fromtensorblock(k)]+=M[k,m]*vj
    end

    r
end

function *{T<:Number}(A::BivariateOperator,b::Vector{T})
    n=size(b,1)
    
    if n>0
        slice(A,:,1:totensorblock(n))*b
    else
        b
    end
end
