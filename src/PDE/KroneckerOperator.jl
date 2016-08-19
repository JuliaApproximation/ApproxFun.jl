export KroneckerOperator



##########
# KroneckerOperator gives the kronecker product of two 1D operators
#########

immutable KroneckerOperator{S,V,DS,RS,DI,RI,T} <: Operator{T}
    ops::Tuple{S,V}
    domainspace::DS
    rangespace::RS
    domaintensorizer::DI
    rangetensorizer::RI
end


KroneckerOperator(A,B,ds::Space,rs::Space,di,ri) =
    KroneckerOperator{typeof(A),typeof(B),typeof(ds),typeof(rs),typeof(di),typeof(ri),
                        promote_type(eltype(A),eltype(B))}((A,B),ds,rs,di,ri)

KroneckerOperator(A,B,ds::Space,rs::Space) = KroneckerOperator(A,B,ds,rs,
                    CachedIterator(tensorizer(ds)),CachedIterator(tensorizer(rs)))
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
    if T == eltype(K)
        K
    else
        ops=Operator{T}(K.ops[1]),Operator{T}(K.ops[2])
        KroneckerOperator{typeof(ops[1]),typeof(ops[2]),typeof(K.domainspace),typeof(K.rangespace),
                            typeof(K.domaintensorizer),typeof(K.rangetensorizer),T}(ops,
              K.domainspace,K.rangespace,
              K.domaintensorizer,K.rangetensorizer)
    end
end


function colstart(A::KroneckerOperator,k::Integer)
    K=block(A.domaintensorizer,k)
    blockstart(A.rangetensorizer,max(1,K-blockbandwidth(A,2)))
end

function colstop(A::KroneckerOperator,k::Integer)
    K=block(A.domaintensorizer,k)
    st=blockstop(A.rangetensorizer,K+blockbandwidth(A,1))
    # zero indicates above dimension
    st==0?size(A,1):min(size(A,1),st)
end

function rowstart(A::KroneckerOperator,k::Integer)
    K=block(rangespace(A),k)
    blockstart(domainspace(A),max(1,K-blockbandwidth(A,1)))
end

function rowstop(A::KroneckerOperator,k::Integer)
    K=block(rangespace(A),k)
    st=blockstop(domainspace(A),K+blockbandwidth(A,2))
    # zero indicates above dimension
    st==0?size(A,2):min(size(A,2),st)
end


bandinds(K::KroneckerOperator) = all(isdiag,K.ops) ? (0,0) : (-∞,∞)

isbandedblockbanded(K::KroneckerOperator) = all(isbanded,K.ops)

blockbandinds(K::KroneckerOperator) =
    bandinds(K.ops[1],1)+bandinds(K.ops[2],1),bandinds(K.ops[1],2)+bandinds(K.ops[2],2)
blockbandinds(K::Operator,k::Integer) = blockbandinds(K)[k]
blockbandwidth(K::Operator,k::Integer) = k==1?-blockbandinds(K,k):blockbandinds(K,k)
subblockbandinds(K::KroneckerOperator,k::Integer) =
    k==1?min(bandinds(K.ops[1],1),-bandinds(K.ops[2],2)):max(bandinds(K.ops[1],2),-bandinds(K.ops[2],1))
subblockbandinds(::Union{ConstantOperator,ZeroOperator},::Integer) = 0


typealias Wrappers Union{ConversionWrapper,MultiplicationWrapper,DerivativeWrapper,LaplacianWrapper,
                       SpaceOperator,ConstantTimesOperator}



isbandedblockbanded(P::Union{PlusOperator,TimesOperator}) = all(isbandedblockbanded,P.ops)



blockbandinds(P::PlusOperator,k::Int) =
    mapreduce(op->blockbandinds(op,k),k==1?min:max,P.ops)
blockbandinds(P::PlusOperator) = blockbandinds(P,1),blockbandinds(P,2)
subblockbandinds(K::PlusOperator,k::Integer) =
    mapreduce(v->subblockbandinds(v,k),k==1?min:max,K.ops)

blockbandinds(P::TimesOperator,k::Int) = mapreduce(op->blockbandinds(op,k),+,P.ops)
subblockbandinds(P::TimesOperator,k::Int) = mapreduce(op->subblockbandinds(op,k),+,P.ops)
blockbandinds(P::TimesOperator) = blockbandinds(P,1),blockbandinds(P,2)

domaintensorizer(P::PlusOperator) = domaintensorizer(P.ops[1])
rangetensorizer(P::PlusOperator) = rangetensorizer(P.ops[1])

domaintensorizer(P::TimesOperator) = domaintensorizer(P.ops[end])
rangetensorizer(P::TimesOperator) = rangetensorizer(P.ops[1])


subblockbandinds(K::Wrappers,k::Integer) = subblockbandinds(K.op,k)
for FUNC in (:blockbandinds,:isbandedblockbanded,:domaintensorizer,:rangetensorizer)
    @eval $FUNC(K::Wrappers) = $FUNC(K.op)
end



function subblockbandindssum(P,k)
    ret=0
    for op in P
        ret+=subblockbandinds(op,k)::Int
    end
    ret
end

subblockbandinds(P::TimesOperator,k)=subblockbandindssum(P.ops,1)

subblockbandinds(K::Operator)=subblockbandinds(K,1),subblockbandinds(K,2)


domainspace(K::KroneckerOperator) = K.domainspace
rangespace(K::KroneckerOperator) = K.rangespace

domaintensorizer(K::KroneckerOperator) = K.domaintensorizer
rangetensorizer(K::KroneckerOperator) = K.rangetensorizer


# we suport 4-indexing with KroneckerOperator
# If A is K x J and B is N x M, then w
# index to match KO=reshape(kron(B,A),N,K,M,J)
# that is
# KO[k,n,j,m] = A[k,j]*B[n,m]
# TODO: arbitrary number of ops

getindex(KO::KroneckerOperator,k::Integer,n::Integer,j::Integer,m::Integer) =
    KO.ops[1][k,j]*KO.ops[2][n,m]

function getindex(KO::KroneckerOperator,kin::Integer,jin::Integer)
    j,m=KO.domaintensorizer[jin]
    k,n=KO.rangetensorizer[kin]
    KO[k,n,j,m]
end

function getindex(KO::KroneckerOperator,k::Integer)
    if size(KO,1) == 1
        KO[1,k]
    elseif size(KO,2) == 1
        KO[k,1]
    else
        throw(ArgumentError("[k] only defined for 1 x ∞ and ∞ x 1 operators"))
    end
end


*(A::KroneckerOperator,B::KroneckerOperator) =
    KroneckerOperator(A.ops[1]*B.ops[1],A.ops[2]*B.ops[2])



## Shorthand


⊗(A,B) = kron(A,B)

Base.kron(A::Operator,B::Operator) = KroneckerOperator(A,B)
Base.kron(A::Operator,B) = KroneckerOperator(A,B)
Base.kron(A,B::Operator) = KroneckerOperator(A,B)
Base.kron{T<:Operator}(A::Vector{T},B::Operator) =
    Operator{promote_type(eltype(T),eltype(B))}[kron(a,B) for a in A]
Base.kron{T<:Operator}(A::Operator,B::Vector{T}) =
    Operator{promote_type(eltype(T),eltype(A))}[kron(A,b) for b in B]
Base.kron{T<:Operator}(A::Vector{T},B::UniformScaling) =
    Operator{promote_type(eltype(T),eltype(B))}[kron(a,1.0B) for a in A]
Base.kron{T<:Operator}(A::UniformScaling,B::Vector{T}) =
    Operator{promote_type(eltype(T),eltype(A))}[kron(1.0A,b) for b in B]


## Conversion




conversion_rule(a::TensorSpace,b::TensorSpace) = conversion_type(a[1],b[1])⊗conversion_type(a[2],b[2])
maxspace(a::TensorSpace,b::TensorSpace) = maxspace(a[1],b[1])⊗maxspace(a[2],b[2])

# TODO: we explicetly state type to avoid type inference bug in 0.4

ConcreteConversion(a::BivariateSpace,b::BivariateSpace) =
    ConcreteConversion{typeof(a),typeof(b),
                        promote_type(eltype(a),eltype(b),real(eltype(domain(a))),real(eltype(domain(b))))}(a,b)

Conversion(a::TensorSpace,b::TensorSpace) = ConversionWrapper(promote_type(eltype(a),eltype(b)),
                KroneckerOperator(Conversion(a[1],b[1]),Conversion(a[2],b[2])))



Multiplication{D,T}(f::Fun{D,T},sp::BivariateSpace) =
    ConcreteMultiplication{D,typeof(sp),T,T}(chop(f,maxabs(f.coefficients)*40*eps(eltype(f))),sp)
function Multiplication{CS<:ConstantSpace,T,V}(f::Fun{TensorSpace{Tuple{CS,V},T,2}},sp::TensorSpace)
    a=Fun(f.coefficients,space(f)[2]) # coefficients are the same
    #Hack to avoid auto-typing bug.  TODO: incorporate basis
    MultiplicationWrapper(eltype(f),f,eye(sp[1])⊗Multiplication(a,sp[2]))
end
function Multiplication{CS<:ConstantSpace,T,V}(f::Fun{TensorSpace{Tuple{V,CS},T,2}},sp::TensorSpace)
    if isempty(f.coefficients)
        a=Fun(zeros(eltype(f),1),space(f)[1])
    else
        a=Fun(f.coefficients,space(f)[1])
    end
    MultiplicationWrapper(eltype(f),f,Multiplication(a,sp[1])⊗eye(sp[2]))
end

Multiplication{D<:UnivariateSpace,SV,TT,T}(f::Fun{D,T},sp::SumSpace{SV,TT,AnyDomain,2})=Multiplication(f⊗1,sp)
Multiplication{D<:UnivariateSpace,T}(f::Fun{D,T},sp::BivariateSpace)=Multiplication(f⊗1,sp)


## PDE Factorization

isfiniterange(::,k) = false
isfiniterange(B::KroneckerOperator,k::Integer) = isfinite(size(B.ops[k],1))
isfiniterange(B::PlusOperator,k::Integer) = isfiniterange(first(B.ops),k)





function findfunctionals(A::Vector,k::Integer)
    T=eltype(eltype(eltype(A)))
    indsBx=find(f->isfiniterange(f,k),A)
    if k==1
        indsBx,Operator{T}[(@assert dekron(Ai,2)==ConstantOperator(Float64,1.0); dekron(Ai,1)) for Ai in A[indsBx]]
    else
        @assert k==2
        indsBx,Operator{T}[(@assert dekron(Ai,1)==ConstantOperator(Float64,1.0); dekron(Ai,2)) for Ai in A[indsBx]]
    end
end




## transpose


Base.transpose(K::KroneckerOperator)=KroneckerOperator(K.ops[2],K.ops[1])

for TYP in (:ConversionWrapper,:MultiplicationWrapper,:DerivativeWrapper,:IntegralWrapper,:LaplacianWrapper),
    FUNC in (:domaintensorizer,:rangetensorizer,:blockbandinds,:subblockbandinds)
    @eval $FUNC(S::$TYP) = $FUNC(S.op)
end


Base.transpose(S::SpaceOperator) =
    SpaceOperator(transpose(S.op),domainspace(S).',rangespace(S).')
Base.transpose(S::ConstantTimesOperator) = sp.c*S.op.'


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
    DerivativeWrapper{typeof(K),typeof(domainspace(K)),Vector{Int},T}(K,order)
end


## Functionals
Evaluation(sp::TensorSpace,x::Vec) = EvaluationWrapper(sp,x,zeros(Int,length(x)),⊗(map(Evaluation,sp.spaces,x)...))
Evaluation(sp::TensorSpace,x::Tuple) = Evaluation(sp,Vec(x...))




### Copy

# finds block lengths for a subrange
function blocklengthrange(rt,kr)
    KR=block(rt,first(kr)):block(rt,last(kr))
    Klengths=Array(Int,length(KR))
    for ν in eachindex(KR)
        Klengths[ν]=blocklength(rt,KR[ν])
    end
    Klengths[1]+=blockstart(rt,KR[1])-kr[1]
    Klengths[end]+=kr[end]-blockstop(rt,KR[end])
    Klengths
end

function Base.convert(::Type{BandedBlockBandedMatrix},S::SubOperator)
    kr,jr=parentindexes(S)
    KO=parent(S)
    l,u=blockbandinds(KO)
    λ,μ=subblockbandinds(KO)

    rt=rangetensorizer(KO)
    dt=domaintensorizer(KO)
    ret=bbbzeros(S)


    for J=1:blocksize(ret,2)
        jshft=jr[1]+blockstart(dt,J)-2
        for K=blockcolrange(ret,J)
            Bs=viewblock(ret,K,J)
            kshft=kr[1]+blockstart(rt,K)-2
            for j=1:size(Bs,2),k=colrange(Bs,j)
                Bs[k,j]=KO[k+kshft,j+jshft]
            end
        end
    end

    ret
end


#TODO: there's a bug here
function Base.convert{KKO<:KroneckerOperator,T}(::Type{BandedBlockBandedMatrix},S::SubOperator{T,KKO})
    kr,jr=parentindexes(S)
    KO=parent(S)
    l,u=blockbandinds(KO)
    λ,μ=subblockbandinds(KO)

    rt=rangetensorizer(KO)
    dt=domaintensorizer(KO)
    ret=bbbzeros(eltype(KO),-l,u,-λ,μ,
                blocklengthrange(rt,kr),
                blocklengthrange(dt,jr))

    A,B=KO.ops
    K=block(rt,kr[end]);J=block(dt,jr[end])
    AA=A[1:K,1:J]
    BB=B[1:K,1:J]



    for J=1:blocksize(ret,2)
        jshft=jr[1]+blockstart(dt,J)-2
        for K=blockcolrange(ret,J)
            Bs=viewblock(ret,K,J)
            kshft=kr[1]+blockstart(rt,K)-2
            for j=1:size(Bs,2),k=colrange(Bs,j)
                κ,ν=subblock2tensor(rt,K,k)
                ξ,μ=subblock2tensor(dt,J,j)
                Bs[k,j]=AA[κ,ξ]*BB[ν,μ]
            end
        end
    end

    ret
end
