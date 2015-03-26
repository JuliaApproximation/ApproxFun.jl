export KroneckerOperator


# gives zero block banded matrix, where the blocks are increasing size
# the bandwidths are assumed to be constant

function blockbandzeros{T}(zer::Function,::Type{T},n,m::Integer,l,u,Bl,Bu)
#    l=Al+Bl;u=Au+Bu
    ret=BandedMatrix(T,n,m,l,u)

    for k=1:n,j=max(1,k-l):min(m,k+u)
#        nl=min(Al,Bu+k-j);nu=min(Au,Bl+j-k)
#        ret[k,j]=zer(eltype(T),k,j,Bl,Bu)
        ret[k,j]=zer(eltype(T),k,j)
    end

    ret
end

blockbandzeros{T}(::Type{T},n,m::Integer,l,u,Bl,Bu)=blockbandzeros(zeros,Matrix{T},n,m,l,u,Bl,Bu)

blockbandzeros{T}(::Type{T},n,m::Colon,Al,Au,Bl,Bu)=blockbandzeros(T,n,n+Au,Al,Au,Bl,Bu)
blockbandzeros{T}(::Type{T},n,m,Alu,Blu)=blockbandzeros(T,n,m,-Alu[1],Alu[2],-Blu[1],Blu[2])
blockbandzeros{T}(zer::Function,::Type{T},n,m,Alu,Blu)=blockbandzeros(zer,T,n,m,-Alu[1],Alu[2],-Blu[1],Blu[2])


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


typealias BivariateOperator{T} BandedOperator{BandedMatrix{T}}


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
KroneckerOperator(A::UniformScaling,B::UniformScaling)=KroneckerOperator(ConstantOperator(A.λ),ConstantOperator(B.λ))
KroneckerOperator(A,B::UniformScaling)=KroneckerOperator(A,ConstantOperator(B.λ))
KroneckerOperator(A::UniformScaling,B)=KroneckerOperator(ConstantOperator(A.λ),B)

for OP in (:promotedomainspace,:promoterangespace)
    @eval $OP(K::KroneckerOperator,ds::TensorSpace)=KroneckerOperator($OP(K.ops[1],ds[1]),
                                                                      $OP(K.ops[2],ds[2]))
end


Base.convert{BO<:Operator}(::Type{BO},K::KroneckerOperator)=KroneckerOperator(convert(BandedOperator{eltype(eltype(BO))},K.ops[1]),
                                                                                convert(BandedOperator{eltype(eltype(BO))},K.ops[2]),
                                                                              K.domainspace,
                                                                              K.rangespace)


bandinds(K::KroneckerOperator)=bandinds(K.ops[1],1)+bandinds(K.ops[2],1),bandinds(K.ops[1],2)+bandinds(K.ops[2],2)
blockbandinds(K::KroneckerOperator,k::Integer)=k==1?min(bandinds(K.ops[1],1),-bandinds(K.ops[2],2)):max(bandinds(K.ops[1],2),-bandinds(K.ops[2],1))

blockbandinds(K::PlusOperator,k::Integer)=mapreduce(v->blockbandinds(v,k),k==1?min:max,K.ops)
blockbandinds(K::ConversionWrapper,k::Integer)=blockbandinds(K.op,k)


blockbandinds{T}(K::BandedOperator{BandedMatrix{T}})=blockbandinds(K,1),blockbandinds(K,2)


domainspace(K::KroneckerOperator)=K.domainspace
rangespace(K::KroneckerOperator)=K.rangespace

function kronaddentries!(A,B,M,kr::Range)
    m=max(size(A,2),size(B,2))
    l=A.l+B.l;u=A.u+B.u

    for k=kr,j=max(1,k-l):k+u
        nl=min(A.l,B.u+k-j);nu=min(A.u,B.l+j-k)
        @inbounds Mkj=M[k,j]
        for κ=1:k,ξ=max(1,κ-nl):min(j,κ+nu)
            #Mkj[κ,ξ]+=A[κ,ξ]*B[k-κ+1,j-ξ+1]
            @inbounds Mkj[κ,ξ]+=A.data[ξ-κ+A.l+1,κ]*B.data[j-k+κ-ξ+B.l+1,k-κ+1]
        end
    end
    M
end

addentries!(K::KroneckerOperator,A,kr::Range)=kronaddentries!(slice(K.ops[1],1:last(kr),:),slice(K.ops[2],1:last(kr),:),A,kr)


bazeros{T}(K::BivariateOperator{T},n::Integer,::Colon)=blockbandzeros(T,n,:,bandinds(K),blockbandinds(K))
bazeros{T}(K::BivariateOperator{T},n::Integer,br::(Int,Int))=blockbandzeros(T,n,:,br,blockbandinds(K))

# function BandedMatrix{T}(K::BivariateOperator{T},kr::UnitRange,::Colon)
#     @assert first(kr)==1
#     BandedMatrix(K,last(kr))
# end


##########
# Multiply a block banded matrix by a vector, where the vector is assumed to
# decompose into blocks
# TODO: Don't assume block banded matrix has i x j blocks
###########

function *{BM<:AbstractArray,V<:Number}(M::BandedMatrix{BM},v::Vector{V})
    n,m=size(M)
    r=zeros(promote_type(eltype(BM),V),div(n*(n+1),2))
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

function *{M,T<:Number}(A::BivariateOperator{M},b::Vector{T})
    n=size(b,1)

    if n>0
        slice(A,:,1:totensorblock(n))*pad(b,fromtensorblock(totensorblock(n))[end])
    else
        b
    end
end


*(A::KroneckerOperator,B::KroneckerOperator)=KroneckerOperator(A.ops[1]*B.ops[1],A.ops[2]*B.ops[2])



## Shorthand


⊗(A,B)=KroneckerOperator(A,B)

Base.kron(A::Operator,B::Operator)=KroneckerOperator(A,B)
Base.kron(A::Operator,B)=KroneckerOperator(A,B)
Base.kron(A,B::Operator)=KroneckerOperator(A,B)

# ⊗{O<:Operator,V<:Operator}(A::Matrix{O},B::Matrix{V})=interlace(A)⊗interlace(B)
# ⊗{O<:Operator}(A::Matrix{O},B::Fun)=interlace(A)⊗B
# ⊗{O<:Operator}(A::Matrix{O},B)=interlace(A)⊗B
# ⊗{O<:Operator}(B::Fun,A::Matrix{O})=A⊗interlace(B)
# ⊗{O<:Operator}(B,A::Matrix{O})=A⊗interlace(B)


## Conversion

conversion_rule(a::TensorSpace,b::TensorSpace)=conversion_type(a[1],b[1])⊗conversion_type(a[2],b[2])
maxspace(a::TensorSpace,b::TensorSpace)=maxspace(a[1],b[1])⊗maxspace(a[2],b[2])

Conversion(a::TensorSpace,b::TensorSpace)=ConversionWrapper(KroneckerOperator(Conversion(a[1],b[1]),Conversion(a[2],b[2])))





## AlmostBandedOperator

function resizedata!{T<:Matrix,M<:BandedOperator,R}(B::AlmostBandedOperator{T,M,R},n::Integer)
    resizedata!(B.fill,n)

    if n > B.datalength
        nbc=B.fill.numbcs

        if n > size(B.data,1)
            newdata=blockbandzeros(eltype(T),n,:,bandinds(B.data),blockbandinds(B.op))

            for k=1:B.datalength,j=max(1,k+bandinds(B.data,1)):k+bandinds(B.data,2)
                newdata.data[k,j]=B.data.data[k,j]
            end
            B.data=newdata
        end

        addentries!(B.op,IndexStride(B.data,nbc,0),B.datalength+1-nbc:n-nbc)
        B.datalength = n
    end

    B
end

function givensmatrix(a::Matrix,b::Matrix)
    q,r=qr([a;b];thin=false)
    q[1:size(a,1),1:size(a,2)],q[1:size(a,1),size(a,1)+1:end],q[size(a,1)+1:end,1:size(a,2)],q[size(a,1)+1:end,size(a,2)+1:end]
end


#TOD: Bcs
function unsafe_getindex{T<:Matrix,R}(B::FillMatrix{T,R},k::Integer,j::Integer)
    zeros(eltype(T),k,j)
end


#TODO: Fix hack override
function pad{T<:Vector}(f::Vector{T},n::Integer)
	if n > length(f)
	   ret=Array(T,n)
	   ret[1:length(f)]=f
	   for j=length(f)+1:n
	       ret[j]=zeros(eltype(T),j)
	   end
       ret
	else
        f[1:n]
	end
end




function backsubstitution!{T<:Vector}(B::AlmostBandedOperator,u::Array{T})
    n=size(u,1)
    b=B.bandinds[end]
    nbc = B.fill.numbcs
    A=B.data

    @assert nbc==0

    for c=1:size(u,2)

        # before we get to filled rows
        for k=n:-1:max(1,n-b)
            @simd for j=k+1:n
                @inbounds u[k,c]-=A.data[j-k+A.l+1,k]*u[j,c]
            end

            @inbounds u[k,c] = A.data[A.l+1,k]\u[k,c]
        end

       #filled rows
        for k=n-b-1:-1:1
            @simd for j=k+1:k+b
                @inbounds u[k,c]-=A.data[j-k+A.l+1,k]*u[j,c]
            end

            @inbounds u[k,c] = A.data[A.l+1,k]\u[k,c]
        end
    end
    u
end



function applygivens!(ca::Matrix,cb,mb,a,B::BandedMatrix,k1::Integer,k2::Integer,jr::Range)
    @simd for j = jr
        @inbounds B1 = B.data[j-k1+B.l+1,k1]    #B[k1,j]
        @inbounds B2 = B.data[j-k2+B.l+1,k2]    #B[k2,j]

        @inbounds B.data[j-k1+B.l+1,k1]=ca*B1
        @inbounds B.data[j-k2+B.l+1,k2]=a*B2
        BLAS.gemm!('N','N',1.0,cb,B2,1.0,B.data[j-k1+B.l+1,k1])
        BLAS.gemm!('N','N',1.0,mb,B1,1.0,B.data[j-k2+B.l+1,k2])
    end

    B
end

function applygivens!(ca::Matrix,cb,mb,a,F::FillMatrix,B::BandedMatrix,k1::Integer,k2::Integer,jr::Range)
    for j = jr
        @inbounds B2 = B.data[j-k2+B.l+1,k2]   #B[k2,j]
        @inbounds B.data[j-k2+B.l+1,k2]=a*B2
    end

    B
end

function applygivens!(ca::Matrix,cb,mb,a,B::Matrix,k1::Integer,k2::Integer)
    for j = 1:size(B,2)
        @inbounds B1 = B[k1,j]
        @inbounds B2 = B[k2,j]

        @inbounds B[k1,j],B[k2,j]= ca*B1 + cb*B2,mb*B1 + a*B2
    end

    B
end
