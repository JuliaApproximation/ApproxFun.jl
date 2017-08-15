export KroneckerOperator



##########
# KroneckerOperator gives the kronecker product of two 1D operators
#########

struct KroneckerOperator{S,V,DS,RS,DI,RI,T} <: Operator{T}
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



function promotedomainspace(K::KroneckerOperator,ds::TensorSpace)
    A=promotedomainspace(K.ops[1],ds.spaces[1])
    B=promotedomainspace(K.ops[2],ds.spaces[2])
    KroneckerOperator(A,B,ds,rangespace(A)⊗rangespace(B))
end

function promoterangespace(K::KroneckerOperator,rs::TensorSpace)
    A=promoterangespace(K.ops[1],rs.spaces[1])
    B=promoterangespace(K.ops[2],rs.spaces[2])
    KroneckerOperator(A,B,domainspace(K),rs)
end


function convert(::Type{Operator{T}},K::KroneckerOperator) where T<:Number
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
    M = blockbandwidth(A,2)
    if isfinite(M)
        blockstart(A.rangetensorizer,max(Block(1),K-M))
    else
        blockstart(A.rangetensorizer,Block(1))
    end
end

function colstop(A::KroneckerOperator,k::Integer)
    k == 0 && return 0
    K=block(A.domaintensorizer,k)
    st=blockstop(A.rangetensorizer,blockcolstop(A,K))
    # zero indicates above dimension
    min(size(A,1),st)
end

function rowstart(A::KroneckerOperator,k::Integer)
    K=block(rangespace(A),k)
    blockstart(domainspace(A),max(Block(1),K-blockbandwidth(A,1)))
end

function rowstop(A::KroneckerOperator,k::Integer)
    K=block(rangespace(A),k)
    st=blockstop(domainspace(A),K+blockbandwidth(A,2))
    # zero indicates above dimension
    st==0?size(A,2):min(size(A,2),st)
end


bandinds(K::KroneckerOperator) = (-∞,∞)

isblockbanded(K::KroneckerOperator) = all(isblockbanded,K.ops)
isbandedblockbanded(K::KroneckerOperator) =
    all(op->isbanded(op) && isinf(size(op,1)) && isinf(size(op,2)),K.ops)
israggedbelow(K::KroneckerOperator) = all(israggedbelow,K.ops)


blockbandinds(K::KroneckerOperator) =
    (blockbandinds(K.ops[1],1)+blockbandinds(K.ops[2],1),
    blockbandinds(K.ops[1],2)+blockbandinds(K.ops[2],2))

# If each block were in turn BlockBandedMatrix, these would
# be the    bandinds
subblock_blockbandinds(K::KroneckerOperator) =
    (min(blockbandinds(K.ops[1],1),-blockbandinds(K.ops[2],2)) ,
           max(blockbandinds(K.ops[1],2),-blockbandinds(K.ops[2],1)))


# If each block were in turn BandedMatrix, these are the bandinds
function subblockbandinds(K::KroneckerOperator)
    if all(hastrivialblocks,domainspace(K).spaces) &&
            all(hastrivialblocks,rangespace(K).spaces)
        subblock_blockbandinds(K)
    else
        dt = domaintensorizer(K).iterator
        rt = rangetensorizer(K).iterator
        # assume block size is repeated and square
        @assert all(b->isa(b,Repeated),dt.blocks)
        @assert rt.blocks == dt.blocks



        sb = subblock_blockbandinds(K)
        # divide by the size of each block
        sb_sz = mapreduce(value,*,dt.blocks)
        # spread by sub block szie
        (sb[1]-1)*sb_sz+1,(sb[2]+1)*sb_sz-1
    end
end

subblockbandinds(K::KroneckerOperator,k::Integer) = subblockbandinds(K)[k]

subblockbandinds(::Union{ConstantOperator,ZeroOperator},::Integer) = 0


const Wrappers = Union{ConversionWrapper,MultiplicationWrapper,DerivativeWrapper,LaplacianWrapper,
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

domaintensorizer(R::Operator) = tensorizer(domainspace(R))
rangetensorizer(R::Operator) = tensorizer(rangespace(R))

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

subblockbandinds(P::TimesOperator,k) = subblockbandindssum(P.ops,1)



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
Base.kron(A::AbstractVector{T},B::Operator) where {T<:Operator} =
    Operator{promote_type(eltype(T),eltype(B))}[kron(a,B) for a in A]
Base.kron(A::Operator,B::AbstractVector{T}) where {T<:Operator} =
    Operator{promote_type(eltype(T),eltype(A))}[kron(A,b) for b in B]
Base.kron(A::AbstractVector{T},B::UniformScaling) where {T<:Operator} =
    Operator{promote_type(eltype(T),eltype(B))}[kron(a,1.0B) for a in A]
Base.kron(A::UniformScaling,B::AbstractVector{T}) where {T<:Operator} =
    Operator{promote_type(eltype(T),eltype(A))}[kron(1.0A,b) for b in B]






## transpose


Base.transpose(K::KroneckerOperator)=KroneckerOperator(K.ops[2],K.ops[1])

for TYP in (:ConversionWrapper,:MultiplicationWrapper,:DerivativeWrapper,:IntegralWrapper,:LaplacianWrapper),
    FUNC in (:domaintensorizer,:rangetensorizer)
    @eval $FUNC(S::$TYP) = $FUNC(S.op)
end


Base.transpose(S::SpaceOperator) =
    SpaceOperator(transpose(S.op),domainspace(S).',rangespace(S).')
Base.transpose(S::ConstantTimesOperator) = sp.c*S.op.'



### Calculus

#TODO: general dimension
function Derivative(S::TensorSpace{SV,DD},order::Vector{Int}) where {SV,DD<:BivariateDomain}
    @assert length(order)==2
    if order[1]==0
        Dy=Derivative(S.spaces[2],order[2])
        K=eye(S.spaces[1])⊗Dy
        T=eltype(Dy)
    elseif order[2]==0
        Dx=Derivative(S.spaces[1],order[1])
        K=Dx⊗eye(S.spaces[2])
        T=eltype(Dx)
    else
        Dx=Derivative(S.spaces[1],order[1])
        Dy=Derivative(S.spaces[2],order[2])
        K=Dx⊗Dy
        T=promote_type(eltype(Dx),eltype(Dy))
    end
    # try to work around type inference
    DerivativeWrapper{typeof(K),typeof(domainspace(K)),Vector{Int},T}(K,order)
end


DefiniteIntegral(S::TensorSpace) = DefiniteIntegralWrapper(mapreduce(DefiniteIntegral,⊗,S.spaces))



### Copy

# finds block lengths for a subrange
blocklengthrange(rt,B::Block) = [blocklength(rt,B)]
blocklengthrange(rt,B::Range{Block}) = blocklength(rt,B)
function blocklengthrange(rt,kr)
    KR=block(rt,first(kr)):block(rt,last(kr))
    Klengths=Vector{Int}(length(KR))
    for ν in eachindex(KR)
        Klengths[ν]=blocklength(rt,KR[ν])
    end
    Klengths[1]+=blockstart(rt,KR[1])-kr[1]
    Klengths[end]+=kr[end]-blockstop(rt,KR[end])
    Klengths
end

function bandedblockbanded_convert!(ret,S::SubOperator,KO,rt,dt)
    pinds = parentindexes(S)
    kr,jr = pinds

    kr1,jr1 = reindex(S,pinds,(1,1))

    Kshft = block(rt,kr1)-1
    Jshft = block(dt,jr1)-1



    for J=Block(1):Block(blocksize(ret,2))
        jshft = (J==Block(1) ? jr1 : blockstart(dt,J+Jshft)) - 1
        for K=blockcolrange(ret,J)
            Bs = view(ret,K,J)
            Bspinds = parentindexes(Bs)
            kshft = (K==Block(1) ? kr1 : blockstart(rt,K+Kshft)) - 1
            for ξ=1:size(Bs,2),κ=colrange(Bs,ξ)
                Bs[κ,ξ] = S[reindex(Bs,Bspinds,(κ,ξ))...]
            end
        end
    end

    ret
end



function default_bandedblockbandedmatrix(S)
    KO = parent(S)
    rt=rangespace(KO)
    dt=domainspace(KO)
    ret=bbbzeros(S)
    bandedblockbanded_convert!(ret,S,parent(S),rt,dt)
end

convert(::Type{BandedBlockBandedMatrix},S::SubOperator) = default_bandedblockbandedmatrix(S)


function convert(::Type{BandedBlockBandedMatrix},
                      S::SubOperator{T,KKO,Tuple{UnitRange{Int},UnitRange{Int}}}) where {KKO<:KroneckerOperator,T}
    kr,jr = parentindexes(S)
    (isempty(kr) || isempty(jr)) && return bbbzeros(S)
    KO = parent(S)

    rt = rangetensorizer(KO)
    dt = domaintensorizer(KO)

    KR = block(rt,kr[1]):block(rt,kr[end])
    JR = block(dt,jr[1]):block(dt,jr[end])

    # use fast block version
    M = BandedBlockBandedMatrix(view(KO,KR,JR))
    k_st = kr[1]-blockstart(rt,KR[1])+1
    j_st = jr[1]-blockstart(dt,JR[1])+1
    ret = M[k_st:k_st+length(kr)-1,j_st:j_st+length(jr)-1]

    # add zero blocks
    prepend!(ret.rows,zeros(Int,KR[1].K-1))
    prepend!(ret.cols,zeros(Int,JR[1].K-1))
    ret.rowblocks[:] += KR[1].K-1
    ret.colblocks[:] += JR[1].K-1
    BSH = (JR[1].K-1) - (KR[1].K-1)
    ret.l -= BSH
    ret.u += BSH

    ret
end


const Trivial2DTensorizer = CachedIterator{Tuple{Int64,Int64},
                                             TrivialTensorizer{2},
                                             Tuple{Tuple{Int64,Int64},Tuple{Int64,Int64},
                                                   Tuple{Int64,Int64},Tuple{Bool,Bool},
                                                   Tuple{Int64,Infinity{Bool}}}}

# This routine is an efficient version of KroneckerOperator for the case of
# tensor product of trivial blocks

function convert(::Type{BandedBlockBandedMatrix},
                      S::SubOperator{T,KroneckerOperator{SS,V,DS,RS,
                                     Trivial2DTensorizer,Trivial2DTensorizer,T},
                                     Tuple{UnitRange{Block},UnitRange{Block}}}) where {SS,V,DS,RS,T}
    KR,JR=parentindexes(S)
    KO=parent(S)

    ret=bbbzeros(S)

    A,B=KO.ops
    AA=A[Block(1):KR[end],Block(1):JR[end]]::BandedMatrix{eltype(S)}
    Al,Au = bandwidths(AA)
    BB=B[Block(1):KR[end],Block(1):JR[end]]::BandedMatrix{eltype(S)}
    Bl,Bu = bandwidths(BB)
    λ,μ = subblockbandwidths(ret)

    for J in Block(1):Block(blocksize(ret,2)), K in blockcolrange(ret,J)
        n,m=KR[K.K].K,JR[J.K].K
        Bs = view(ret,K,J)
        l = min(Al,Bu+n-m,λ)
        u = min(Au,Bl+m-n,μ)
        @inbounds for j=1:m, k=max(1,j-u):min(n,j+l)
            a = inbands_getindex(AA,k,j)
            b = inbands_getindex(BB,n-k+1,m-j+1)
            c = a*b
            inbands_setindex!(Bs,c,k,j)
        end
    end
    ret
end

## TensorSpace operators


## Conversion




conversion_rule(a::TensorSpace,b::TensorSpace) = conversion_type(a.spaces[1],b.spaces[1])⊗conversion_type(a.spaces[2],b.spaces[2])
maxspace(a::TensorSpace,b::TensorSpace) = maxspace(a.spaces[1],b.spaces[1])⊗maxspace(a.spaces[2],b.spaces[2])

# TODO: we explicetly state type to avoid type inference bug in 0.4

ConcreteConversion(a::BivariateSpace,b::BivariateSpace) =
    ConcreteConversion{typeof(a),typeof(b),
                        promote_type(prectype(a),prectype(b))}(a,b)

Conversion(a::TensorSpace,b::TensorSpace) = ConversionWrapper(promote_type(prectype(a),prectype(b)),
                KroneckerOperator(Conversion(a.spaces[1],b.spaces[1]),Conversion(a.spaces[2],b.spaces[2])))



function Multiplication(f::Fun{TS},S::TensorSpace) where {TS<:TensorSpace}
    lr=LowRankFun(f)
    ops=map(kron,map(a->Multiplication(a,S.spaces[1]),lr.A),map(a->Multiplication(a,S.spaces[2]),lr.B))
    MultiplicationWrapper(f,+(ops...))
end

## Functionals
Evaluation(sp::TensorSpace,x::Vec) = EvaluationWrapper(sp,x,zeros(Int,length(x)),⊗(map(Evaluation,sp.spaces,x)...))
Evaluation(sp::TensorSpace,x::Tuple) = Evaluation(sp,Vec(x...))



# it's faster to build the operators to the last b
function A_mul_B_coefficients(A::SubOperator{T,KKO,Tuple{UnitRange{Int},UnitRange{Int}}},b) where {T,KKO<:KroneckerOperator}
    P = parent(A)
    kr,jr = parentindexes(A)
    dt,rt = domaintensorizer(P),rangetensorizer(P)
    KR,JR = Block(1):block(rt,kr[end]),Block(1):block(dt,jr[end])
    M = P[KR,JR]
    view(M,kr,jr)*b
end
