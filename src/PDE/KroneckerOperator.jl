export KroneckerOperator



##########
# KroneckerOperator gives the kronecker product of two 1D operators
#########

immutable KroneckerOperator{S,V,DS,RS,T<:Number} <: Operator{T}
    ops::Tuple{S,V}
    domainspace::DS
    rangespace::RS
end


KroneckerOperator(A,B,ds,rs) =
    KroneckerOperator{typeof(A),typeof(B),typeof(ds),typeof(rs),
                        promote_type(eltype(A),eltype(B))}((A,B),ds,rs)
function KroneckerOperator(A,B)
    ds=domainspace(A)⊗domainspace(B)
    rs=rangespace(A)⊗rangespace(B)
    KroneckerOperator(A,B,ds,rs)
end
KroneckerOperator(A::UniformScaling,B::UniformScaling) =
    KroneckerOperator(ConstantOperator(A.λ),ConstantOperator(B.λ))
KroneckerOperator(A,B::UniformScaling) = KroneckerOperator(A,ConstantOperator(B.λ))
KroneckerOperator(A::UniformScaling,B) = KroneckerOperator(ConstantOperator(A.λ),B)
KroneckerOperator(A::Fun,B::Fun) = KroneckerOperator(Multiplication(A),Multiplication(B))
KroneckerOperator(A::UniformScaling,B::Fun) = KroneckerOperator(ConstantOperator(A.λ),Multiplication(B))
KroneckerOperator(A::Fun,B::UniformScaling) = KroneckerOperator(Multiplication(A),ConstantOperator(B.λ))
KroneckerOperator(A,B::Fun) = KroneckerOperator(A,Multiplication(B))
KroneckerOperator(A::Fun,B) = KroneckerOperator(Multiplication(A),B)



Base.eye{T,D}(S::Space{T,D,2}) = KroneckerOperator(IdentityOperator(),IdentityOperator(),S,S)


function promotedomainspace(K::KroneckerOperator,ds::TensorSpace)
    A=promotedomainspace(K.ops[1],ds[1])
    B=promotedomainspace(K.ops[2],ds[2])
    KroneckerOperator(A,B,ds,rangespace(A)⊗rangespace(B))
end

function promoterangespace(K::KroneckerOperator,rs::TensorSpace)
    A=promoterangespace(K.ops[1],rs[1])
    B=promoterangespace(K.ops[2],rs[2])
    KroneckerOperator(A,B,domainspace(K),rs)
end


function Base.convert{T<:Number}(::Type{Operator{T}},K::KroneckerOperator)
    if BandedMatrix{T} == eltype(K)
        K
    else
        ops=convert(Operator{T},K.ops[1]),convert(Operator{T},K.ops[2])
        KroneckerOperator{typeof(ops[1]),typeof(ops[2]),typeof(K.domainspace),typeof(K.rangespace),T}(ops,
              K.domainspace,
              K.rangespace)
    end
end


blockbandinds(K::KroneckerOperator) =
    bandinds(K.ops[1],1)+bandinds(K.ops[2],1),bandinds(K.ops[1],2)+bandinds(K.ops[2],2)
subblockbandinds(K::KroneckerOperator,k::Integer) =
    k==1?min(bandinds(K.ops[1],1),-bandinds(K.ops[2],2)):max(bandinds(K.ops[1],2),-bandinds(K.ops[2],1))
subblockbandinds(::Union{ConstantOperator,ZeroOperator},::Integer)=0
subblockbandinds(K::Union{ConversionWrapper,MultiplicationWrapper,DerivativeWrapper,LaplacianWrapper,
                       SpaceOperator,ConstantTimesOperator},k::Integer)=subblockbandinds(K.op,k)



subblockbandinds(K::PlusOperator,k::Integer) =
    mapreduce(v->subblockbandinds(v,k),k==1?min:max,K.ops)

Base.copy{T<:BandedMatrix}(P::SubBandedMatrix{T,PlusOperator{T}}) =
        default_copy(P)


function subblockbandindssum(P,k)
    ret=0
    for op in P
        ret+=subblockbandinds(op,k)::Int
    end
    ret
end

subblockbandinds(P::TimesOperator,k)=subblockbandindssum(P.ops,1)

subblockbandinds{BT<:BandedMatrix}(K::Operator{BT})=subblockbandinds(K,1),subblockbandinds(K,2)


for OP in (:domainspace,:rangespace)
    @eval $OP{BT<:BandedMatrix}(K::Operator{BT},k::Integer)=$OP(K)[k]
end
domainspace(K::KroneckerOperator) = K.domainspace
rangespace(K::KroneckerOperator) = K.rangespace

function getindex(K::KroneckerOperator,kin::Integer,jin::Integer)
    k=totensorblock(kin)
    j=totensorblock(jin)
    T=eltype(K)
    
    if bandinds(K,1) ≤ j-k ≤ bandinds(K,2)
        A=K.ops[1][1:k,1:j]
        B=K.ops[2][1:k,1:j]
        nl=max(0,min(A.l,B.u+k-j));nu=max(0,min(A.u,B.l+j-k))
        ret=BandedMatrix(T,k,j,nl,nu)
        for (κ,ξ) in eachbandedindex(ret)
            ret[κ,ξ]=A[κ,ξ]*B[k-κ+1,j-ξ+1]
        end
        ret
    else
        bzeros(T,k,j,0,0)
    end
end


bzeros{BT<:BandedMatrix}(K::Operator{BT},
                          n::Integer,
                          ::Colon) = blockbandzeros(eltype(BT),n,:,bandinds(K),subblockbandinds(K))
bzeros{BT<:BandedMatrix}(K::Operator{BT},n::Integer,
                          br::Tuple{Int,Int}) = error("Fix call signature") #blockbandzeros(eltype(BT),n,:,br,subblockbandinds(K))


isbzeros{BT<:BandedMatrix}(K::Operator{BT},
                          rws::Range,
                          ::Colon) = isblockbandzeros(eltype(BT),rws,:,bandinds(K),subblockbandinds(K))

bzeros{BT<:BandedMatrix}(K::Operator{BT},rws::Range,
                          br::Tuple{Int,Int}) = error("Fix call signature")   # isblockbandzeros(eltype(BT),n,:,br,subblockbandinds(K))


# function BandedMatrix{T}(K::BivariateOperator{T},kr::UnitRange,::Colon)
#     @assert first(kr)==1
#     BandedMatrix(K,last(kr))
# end





##########
# Multiply a block banded matrix by a vector, where the vector is assumed to
# decompose into blocks
# TODO: Don't assume block banded matrix has i x j blocks
###########

function *{BM<:AbstractArray,TT<:Number}(M::BandedMatrix{BM},v::Vector{TT})
    n,m=size(M)
    r=zeros(promote_type(eltype(BM),eltype(v)),div(n*(n+1),2))
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

function *{BM<:BandedMatrix}(A::Operator{BM},b::Vector)
    n=size(b,1)

    ret=if n>0
        BandedMatrix(A,:,1:totensorblock(n))*pad(b,fromtensorblock(totensorblock(n))[end])
    else
        b
    end

    Fun(ret,rangespace(A))
end


*(A::KroneckerOperator,B::KroneckerOperator)=KroneckerOperator(A.ops[1]*B.ops[1],A.ops[2]*B.ops[2])



## Shorthand


⊗(A,B) = kron(A,B)

Base.kron(A::Operator,B::Operator) = KroneckerOperator(A,B)
Base.kron(A::Operator,B) = KroneckerOperator(A,B)
Base.kron(A,B::Operator) = KroneckerOperator(A,B)
Base.kron{T<:Operator}(A::Vector{T},B::Operator) =
    Operator{BandedMatrix{promote_type(eltype(T),eltype(B))}}[kron(a,B) for a in A]
Base.kron{T<:Operator}(A::Operator,B::Vector{T}) =
    Operator{BandedMatrix{promote_type(eltype(T),eltype(A))}}[kron(A,b) for b in B]
Base.kron{T<:Operator}(A::Vector{T},B::UniformScaling) =
    Operator{BandedMatrix{promote_type(eltype(T),eltype(B))}}[kron(a,1.0B) for a in A]
Base.kron{T<:Operator}(A::UniformScaling,B::Vector{T}) =
    Operator{BandedMatrix{promote_type(eltype(T),eltype(A))}}[kron(1.0A,b) for b in B]


## Conversion




conversion_rule(a::TensorSpace,b::TensorSpace) = conversion_type(a[1],b[1])⊗conversion_type(a[2],b[2])
conversion_rule(b::TensorSpace{AnySpace,AnySpace},a::TensorSpace) = a
conversion_rule(b::TensorSpace{AnySpace,AnySpace},a::Space) = a
maxspace(a::TensorSpace,b::TensorSpace) = maxspace(a[1],b[1])⊗maxspace(a[2],b[2])

# TODO: we explicetly state type to avoid type inference bug in 0.4

ConcreteConversion(a::BivariateSpace,b::BivariateSpace) =
    ConcreteConversion{typeof(a),typeof(b),
                        BandedMatrix{promote_type(eltype(a),eltype(b),real(eltype(domain(a))),real(eltype(domain(b))))}}(a,b)

Conversion(a::TensorSpace,b::TensorSpace) = ConversionWrapper(BandedMatrix{promote_type(eltype(a),eltype(b))},
                KroneckerOperator(Conversion(a[1],b[1]),Conversion(a[2],b[2])))


function Conversion(a::BivariateSpace,b::BivariateSpace)
    if a==b
        error("Don't call conversion to itself")
    elseif conversion_type(a,b)==NoSpace()
        sp=canonicalspace(a)
        if typeof(sp) == typeof(a)
            error("Implement Conversion from " * string(typeof(sp)) * " to " * string(typeof(b)))
        elseif typeof(sp) == typeof(b)
            error("Implement Conversion from " * string(typeof(a)) * " to " * string(typeof(sp)))
        else
            Conversion(a,sp,b)
        end
    else
        error("Implement Conversion from " * string(typeof(a)) * " to " * string(typeof(b)))
    end
end



Multiplication{D,T}(f::Fun{D,T},sp::BivariateSpace) =
    ConcreteMultiplication{D,typeof(sp),T,BandedMatrix{T}}(chop(f,maxabs(f.coefficients)*40*eps(eltype(f))),sp)
function Multiplication{CS<:ConstantSpace,T,V}(f::Fun{TensorSpace{Tuple{CS,V},T,2}},sp::TensorSpace)
    a=Fun(vec(totensor(f.coefficients)[1,:]),space(f)[2])
    #Hack to avoid auto-typing bug.  TODO: incorporate basis
    MultiplicationWrapper(BandedMatrix{eltype(f)},f,eye(sp[1])⊗Multiplication(a,sp[2]))
end
function Multiplication{CS<:ConstantSpace,T,V}(f::Fun{TensorSpace{Tuple{V,CS},T,2}},sp::TensorSpace)
    if isempty(f.coefficients)
        a=Fun(zeros(eltype(f),1),space(f)[1])
    else
        a=Fun(totensor(f.coefficients)[:,1],space(f)[1])
    end
    MultiplicationWrapper(BandedMatrix{eltype(f)},f,Multiplication(a,sp[1])⊗eye(sp[2]))
end

Multiplication{D<:UnivariateSpace,SV,TT,T}(f::Fun{D,T},sp::SumSpace{SV,TT,AnyDomain,2})=Multiplication(f⊗1,sp)
Multiplication{D<:UnivariateSpace,T}(f::Fun{D,T},sp::BivariateSpace)=Multiplication(f⊗1,sp)



# from algebra
function promotedomainspace{T,T2}(P::PlusOperator{T},sp::Space,cursp::TensorSpace{AnySpace,AnySpace,T2})
    if sp==cursp
        P
    else
        promoteplus(Operator{T}[promotedomainspace(op,sp) for op in P.ops])
    end
end

function promotedomainspace{T}(P::TimesOperator,sp::Space,cursp::TensorSpace{AnySpace,AnySpace,T})
    if sp==cursp
        P
    elseif length(P.ops)==2
        P.ops[1]*promotedomainspace(P.ops[end],sp)
    else
        TimesOperator(P.ops[1:end-1])*promotedomainspace(P.ops[end],sp)
    end
end

for op in (:promotedomainspace,:promoterangespace)
    @eval $op(P::Operator,sp::Space,::TensorSpace{AnySpace,AnySpace})=SpaceOperator(P,sp)
end


## PDE Factorization

isafunctional(::,k) = false
isafunctional(B::KroneckerOperator,k::Integer) = isafunctional(B.ops[k])
isafunctional(B::PlusOperator,k::Integer) = isafunctional(first(B.ops),k)





function findfunctionals(A::Vector,k::Integer)
    T=eltype(eltype(eltype(A)))
    indsBx=find(f->isafunctional(f,k),A)
    if k==1
        indsBx,Operator{T}[(@assert dekron(Ai,2)==ConstantOperator(Float64,1.0); dekron(Ai,1)) for Ai in A[indsBx]]
    else
        @assert k==2
        indsBx,Operator{T}[(@assert dekron(Ai,1)==ConstantOperator(Float64,1.0); dekron(Ai,2)) for Ai in A[indsBx]]
    end
end




## transpose

function transpose{T}(A::PlusOperator{BandedMatrix{T}})
    @assert all(map(iskronop,A.ops))
    PlusOperator(Operator{BandedMatrix{eltype(eltype(A))}}[op.' for op in A.ops])
end


Base.transpose(K::KroneckerOperator)=KroneckerOperator(K.ops[2],K.ops[1])

for TYP in (:ConversionWrapper,:MultiplicationWrapper,:DerivativeWrapper,:IntegralWrapper,:LaplacianWrapper)
    @eval Base.transpose(S::$TYP)=$TYP(transpose(S.op))
end


Base.transpose(S::SpaceOperator) =
    SpaceOperator(transpose(S.op),domainspace(S).',rangespace(S).')
Base.transpose(S::ConstantTimesOperator) = sp.c*S.op.'
Base.transpose{V,T<:AbstractArray}(C::ConstantOperator{V,T},k) = C
Base.transpose{V,T<:AbstractArray}(C::ZeroOperator{V,T},k) = C



### Calculus

#TODO: general dimension
function Derivative{SV,TT}(S::TensorSpace{SV,TT,2},order::Vector{Int})
    @assert length(order)==2
    if order[1]==0
        Dy=Derivative(S[2],order[2])
        K=eye(S[1])⊗Dy
        T=eltype(Dy)
    elseif order[2]==0
        Dx=Derivative(S[1],order[1])
        K=Dx⊗eye(S[2])
        T=eltype(Dx)
    else
        Dx=Derivative(S[1],order[1])
        Dy=Derivative(S[2],order[2])
        K=Dx⊗Dy
        T=promote_type(eltype(Dx),eltype(Dy))
    end
    # try to work around type inference
    DerivativeWrapper{typeof(K),typeof(domainspace(K)),Vector{Int},BandedMatrix{T}}(K,order)
end
