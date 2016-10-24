export ⊕,depiece,pieces,PiecewiseSpace



## BlockInterlacer
# interlaces coefficients by blocsk
#

immutable BlockInterlacer{DMS<:Tuple}
    blocks::DMS
end

BlockInterlacer(v::Vector) = BlockInterlacer(tuple(v...))

Base.eltype(it::BlockInterlacer) = Tuple{Int,Int}

dimensions(b::BlockInterlacer) = map(length,b.blocks)
Base.length(b::BlockInterlacer) = mapreduce(length,+,b.blocks)

# the state is always (whichblock,curblock,cursubblock,curcoefficients)
Base.start(it::BlockInterlacer) = (1,1,map(start,it.blocks),ntuple(zero,length(it.blocks)))

function Base.next(it::BlockInterlacer,st)
    N,k,blkst,lngs = st

    if N>length(it.blocks)
        # increment to next block
        blkst = map((blit,blst)->done(blit,blst)?blst:next(blit,blst)[2],it.blocks,blkst)
        return next(it,(1,1,blkst,lngs))
    end

    if done(it.blocks[N],blkst[N])
        # increment to next N
        return next(it,(N+1,1,blkst,lngs))
    end

    B,nxtb = next(it.blocks[N],blkst[N])  # B is block size

    if k > B
        #increment to next N
        return next(it,(N+1,1,blkst,lngs))
    end


    lngs = tuple(lngs[1:N-1]...,lngs[N]+1,lngs[N+1:end]...)
    return (N,lngs[N]),(N,k+1,blkst,lngs)
end


# are all Ints, so finite dimensional
function Base.done(it::BlockInterlacer,st)
    for k=1:length(it.blocks)
        if st[end][k] ≤ length(it.blocks[k])
            return false
        end
    end
    return true
end


#
# immutable RestrictInterlacer{III}
#     interlacer::III
#     sub::Int
# end
#
# Base.eltype{III}(R::Type{RestrictInterlacer{III}) = Int
# Base.eltype(R::RestrictInterlacer) = Int
# Base.start(R::RestrictInterlacer) = (1,start(R.interlacer))
#
# function Base.next(R::RestrictInterlacer,st)
#     a,nx = next(R.interlacer,st[2])
#     if a[1] == R.K
#         st[1],(st[1]+1,nx)
#     else
#         next(st[
#
# end


function findind(it::BlockInterlacer,K::Int,kr::Range)

end



## SumSpace encodes a space that can be decoupled as f(x) = a(x) + b(x) where a is in S and b is in V


abstract DirectSumSpace{SV,T,DD,d} <: Space{T,DD,d}


dimension(sp::DirectSumSpace) = mapreduce(dimension,+,sp.spaces)

spaces(s::Space) = (s,)
spaces(sp::DirectSumSpace) = sp.spaces
space(s::Space,k...) = spaces(s)[k...]
space(f::Fun,k...) = space(space(f),k...)

BlockInterlacer(sp::DirectSumSpace) = BlockInterlacer(map(blocklengths,sp.spaces))
interlacer(sp::DirectSumSpace) = BlockInterlacer(sp)
interlacer(sp::Space) = BlockInterlacer(tuple(blocklengths(sp)))
cache(Q::BlockInterlacer) = CachedIterator(Q)

function blocklengths(sp::DirectSumSpace)
    bl=map(blocklengths,sp)
    N=mapreduce(length,max,bl)
    mapreduce(b->pad(b,N),+,bl)
end
block(sp::DirectSumSpace,k::Int) = findfirst(x->x≥k,cumsum(blocklengths(sp)))



for TYP in (:SumSpace,:TupleSpace)
    @eval begin
        immutable $TYP{SV,T,DD,d} <: DirectSumSpace{SV,T,DD,d}
            spaces::SV
            $TYP(dom::Domain) = new(tuple(map(typ->typ(dom),SV.parameters)...))
            $TYP(sp::Tuple) = new(sp)
        end

        $TYP(sp::Tuple) = $TYP{typeof(sp),mapreduce(basistype,promote_type,sp),
                             typeof(domain(first(sp))),domaindimension(first(sp))}(sp)
    end
end

immutable PiecewiseSpace{SV,T,DD<:UnionDomain,d} <: DirectSumSpace{SV,T,DD,d}
    spaces::SV
    PiecewiseSpace(dom::AnyDomain)=new(tuple(map(typ->typ(dom),SV.parameters)...))
    PiecewiseSpace(dom::UnionDomain)=new(tuple(map((typ,dom)->typ(dom),SV.parameters,dom.domains)...))
    PiecewiseSpace(sp::Tuple)=new(sp)
end

function PiecewiseSpace(spin::Tuple)
    sp=tuple(union(spin)...)  # remove duplicates

    PiecewiseSpace{typeof(sp),mapreduce(basistype,promote_type,sp),
                               typeof(UnionDomain(map(domain,sp))),domaindimension(first(sp))}(sp)
end



for TYP in (:SumSpace,:TupleSpace,:PiecewiseSpace)
    @eval begin
        $TYP(A::$TYP,B::$TYP)=$TYP(tuple(A.spaces...,B.spaces...))

        $TYP(A::Space,B::$TYP)=$TYP(tuple(A,B.spaces...))
        $TYP(A::$TYP,B::Space)=$TYP(tuple(A.spaces...,B))
        $TYP(A::Space...)=$TYP(A)
        $TYP(sp::Array)=$TYP(tuple(sp...))

        canonicalspace(A::$TYP)=$TYP(sort([A.spaces...]))
    end
end


for TYP in (:SumSpace,:TupleSpace)
    @eval setdomain(A::$TYP,d::Domain)=$TYP(map(sp->setdomain(sp,d),A.spaces))
end

setdomain(A::PiecewiseSpace,d::UnionDomain)=PiecewiseSpace(map((sp,dd)->setdomain(sp,dd),A.spaces,d.domains))



function spacescompatible{S<:DirectSumSpace}(A::S,B::S)
    if length(A) != length(B)
        false
    else
        ret=true
        for k=1:length(A)
            if !spacescompatible(A[k],B[k])
                return false
            end
        end
        true
    end
end

function Base.promote_rule{SV,B,DD,d,V,T<:Number}(::Type{Fun{SumSpace{SV,B,DD,d},V}},::Type{T})
    for k=1:length(SV.parameters)
        pt=promote_type(Fun{SV.parameters[k],V},T)
        if pt != Fun
            return Fun{SumSpace{Tuple{SV.parameters[1:k-1]...,pt.parameters[1],SV.parameters[k+1:end]...},
                       B,DD,d},promote_type(V,T)}
        end
    end
    Fun
end

Base.promote_rule{SV,B,DD,d,T<:Number}(::Type{Fun{SumSpace{SV,B,DD,d}}},::Type{T})=promote_rule(Fun{SumSpace{SV,B,DD,d},Float64},T)

function Base.promote_rule{SV,B,DD,d,V,T<:Number}(::Type{Fun{PiecewiseSpace{SV,B,DD,d},V}},::Type{T})
    # if any doesn't support promoting, just leave unpromoted

    newfsp=map(s->promote_type(Fun{s,V},T),SV.parameters)
    if any(s->s==Fun,newfsp)
        Fun
    else
        newsp=map(s->s.parameters[1],newfsp)
        Fun{PiecewiseSpace{Tuple{newsp...},B,DD,d},promote_type(V,T)}
    end
end

Base.promote_rule{SV,B,DD,d,T<:Number}(::Type{Fun{PiecewiseSpace{SV,B,DD,d}}},::Type{T})=promote_rule(Fun{PiecewiseSpace{SV,B,DD,d},Float64},T)



for OP in (:(Base.length),:(Base.start),:(Base.endof))
    @eval $OP(S::DirectSumSpace)=$OP(spaces(S))
end

for OP in (:(Base.getindex),:(Base.next),:(Base.done))
    @eval $OP(S::DirectSumSpace,k)=$OP(spaces(S),k)
end


#support tuple set
for OP in (:(Base.start),:(Base.done),:(Base.endof))
    @eval $OP{SS<:DirectSumSpace}(f::Fun{SS},k...)=$OP(space(f),k...)
end

Base.next{SS<:DirectSumSpace}(f::Fun{SS},k)=f[k],k+1


for TYP in (:SumSpace,:TupleSpace)
    @eval domain(A::$TYP)=domain(A.spaces[end])      # TODO: this assumes all spaces have the same domain
                                                     #        we use end to avoid ConstantSpace
end

Space(d::UnionDomain)=PiecewiseSpace(map(Space,d.domains))
domain(S::PiecewiseSpace)=UnionDomain(map(domain,S.spaces))



function spacescompatible(A::Tuple,B::Tuple)
    if length(A) != length(B)
        return false
    end
    #assumes domain doesn't impact sorting
    asort=sort([A...]);bsort=sort([B...])
    for k=1:length(asort)
        if !spacescompatible(asort[k],bsort[k])
            return false
        end
    end

    return true
end


# avoids default ConstantSpace
function union_rule(B::ConstantSpace,A::SumSpace)
    if !domainscompatible(A,B)
        NoSpace()
    else
        for sp in A.spaces
            if isconvertible(B,sp)
                return A
            end
        end
        SumSpace(A,B)
    end
end

function union_rule(A::SumSpace,B::Space)
    if !domainscompatible(A,B)
        NoSpace()
    else
        for sp in A.spaces
            if isconvertible(B,sp)
                return A
            end
        end
        SumSpace(A,B)
    end
end



## evaluate


for OP in (:(Base.last),:(Base.first))
    @eval begin
        $OP{SS<:SumSpace}(f::Fun{SS})=mapreduce($OP,+,vec(f))
        $OP{SS<:PiecewiseSpace}(f::Fun{SS})=$OP($OP(vec(f)))
        $OP{SS<:TupleSpace}(f::Fun{SS})=$OP(vec(f))
    end
end

evaluate(f::AbstractVector,S::SumSpace,x)=mapreduce(vf->evaluate(vf,x),+,vec(Fun(f,S)))


function evaluate(f::AbstractVector,S::PiecewiseSpace,x)
    d=domain(S)
    g=Fun(f,S)

#    ret=zero(promote_type(eltype(f),eltype(S)))
    for k=1:numpieces(d)
        if in(x,d[k])
            return g[k](x)
        end
    end
    return 0*first(g)
end

function evaluate(v::AbstractVector,S::PiecewiseSpace,x::AbstractVector)
    f=Fun(v,S)
    [f(xk) for xk in x]
end

function evaluate(v::AbstractVector,S::TupleSpace,x...)
    f=Fun(v,S)
    eltype(f)[f[k](x...) for k=1:length(f.space)]
end


## calculus
for TYP in (:SumSpace,:TupleSpace,:PiecewiseSpace)
    for OP in (:differentiate,:integrate)
        @eval function $OP{D<:$TYP,T}(f::Fun{D,T})
            fs=map($OP,vec(f))
            sp=$TYP(map(space,fs))
            Fun(interlace(fs,sp),sp)
        end
    end
    for OP in (:(Base.real),:(Base.imag),:(Base.conj))
        @eval begin
            $OP{SV,DD,d}(f::Fun{$TYP{SV,RealBasis,DD,d}}) = Fun($OP(f.coefficients),f.space)
            function $OP{SV,T,DD,d}(f::Fun{$TYP{SV,T,DD,d}})
                fs=map($OP,vec(f))
                sp=$TYP(map(space,fs))
                Fun(interlace(fs,sp),sp)
            end
        end
    end
end


for TYP in (:SumSpace,:TupleSpace)
    @eval function Base.cumsum{D<:$TYP,T}(f::Fun{D,T})
        fs=map(cumsum,vec(f))
        sp=$TYP(map(space,fs))
        Fun(interlace(fs,sp),sp)
    end
end


for TYP in (:SumSpace,:PiecewiseSpace), OP in (:(Base.sum),:linesum)
    @eval $OP{V<:$TYP}(f::Fun{V})=mapreduce($OP,+,vec(f))
end

function Base.cumsum{V<:PiecewiseSpace}(f::Fun{V})
    vf=pieces(f)
    r=zero(eltype(f))
    for k=1:length(vf)
        vf[k]=cumsum(vf[k]) + r
        r=last(vf[k])
    end
    depiece(vf)
end

Base.cumsum{V<:PiecewiseSpace}(f::Fun{V},d::Domain)=mapreduce(g->cumsum(g,d),+,pieces(f))



bilinearform{S<:PiecewiseSpace,V<:PiecewiseSpace}(f::Fun{S},g::Fun{V}) = sum(map(bilinearform,pieces(f),pieces(g)))
linebilinearform{S<:PiecewiseSpace,V<:PiecewiseSpace}(f::Fun{S},g::Fun{V}) = sum(map(linebilinearform,pieces(f),pieces(g)))

# assume first domain has 1 as a basis element



function Base.ones{T<:Number}(::Type{T},S::SumSpace)
    @assert length(S.spaces)==2
    if isconvertible(ConstantSpace(),S.spaces[1])
        ones(T,S[1])⊕zeros(T,S[2])
    else
        zeros(T,S[1])⊕ones(T,S[2])
    end
end

Base.ones(S::SumSpace)=ones(Float64,S)

Base.ones{T<:Number,SS,V}(::Type{T},S::PiecewiseSpace{SS,V}) = depiece(map(ones,spaces(S)))
Base.ones(S::PiecewiseSpace) = ones(Float64,S)


identity_fun(S::PiecewiseSpace) = depiece(map(identity_fun,S.spaces))

# vec

function Base.getindex{DSS<:DirectSumSpace}(f::Fun{DSS},k::Integer)
    it=interlacer(space(f))
    N=length(f.coefficients)
    d=dimension(space(f,k))

    # preallocate: we know we have at most N coefficients
    ret=Array(eltype(f),N)
    j=1  # current coefficient
    p=0  # current length
    for (n,m) in it
        if j > N
            break
        end
        if n==k
            ret[m]=f.coefficients[j]
            p+=1
            if m ≥ d
                # if we've reached the dimension, we are done
                break
            end
        end
        j+=1
    end
    resize!(ret,p)  # through out extra coefficients
    Fun(ret,space(f,k))
end


# interlace coefficients according to iterator
function interlace{T,V<:AbstractVector}(::Type{T},v::AbstractVector{V},it::BlockInterlacer)
    ret=Vector{T}()
    N=mapreduce(length,max,v)
    cnts = map(length,v)

    for (n,m) in it
        if maximum(cnts) == 0
            break
        end

        if m ≤ length(v[n])
            push!(ret,v[n][m])
            cnts[n] -= 1
        else
            push!(ret,zero(T))
        end
    end
    ret
end

interlace{V<:AbstractVector}(v::AbstractVector{V},it::BlockInterlacer) =
    interlace(mapreduce(eltype,promote_type,v),v,it)

interlace{V<:AbstractVector}(v::AbstractVector{V},sp::DirectSumSpace) =
    interlace(v,interlacer(sp))

interlace{F<:Fun}(v::AbstractVector{F},sp::DirectSumSpace) =
    interlace(map(coefficients,v),sp)


function interlace(v::Union{Tuple,Vector{Any}},sp::DirectSumSpace)
    if all(vk->isa(vk,Fun),v)
        V=Array(Vector{mapreduce(eltype,promote_type,v)},length(v))
        for k=1:length(v)
            V[k]=coefficients(v[k])
        end
    else
        V=Array(Vector{mapreduce(eltype,promote_type,v)},length(v))
        for k=1:length(v)
            V[k]=v[k]
        end
    end
    interlace(V,sp)
end


Base.vec(S::DirectSumSpace) = S.spaces
Base.vec{S<:DirectSumSpace}(f::Fun{S}) = Fun[f[j] for j=1:length(f.space)]

pieces{S<:PiecewiseSpace}(f::Fun{S}) = vec(f)

for (Dep,Sp) in ((:depiece,:PiecewiseSpace),(:detuple,:TupleSpace))
    @eval begin
        function $Dep{F<:Fun}(v::Vector{F})
            spaces=map(space,v)
            sp=$Sp(spaces)
            Fun(interlace(v,sp),sp)
        end
        function $Dep(v::Tuple)
            spaces=map(space,v)
            sp=$Sp(spaces)
            Fun(interlace(v,sp),sp)
        end

        $Dep(v::Vector{Any})=$Dep(tuple(v...))
    end
end

interlace{FF<:Fun}(f::AbstractVector{FF}) = vcat(f...)

# convert a vector to a Fun with TupleSpace
Fun(v::Vector{Any},sp::TupleSpace) = detuple(map(Fun,v,sp.spaces))
Fun{F<:Fun}(v::Vector{F},sp::TupleSpace) = detuple(map(Fun,v,sp.spaces))



## transforms


function points(d::PiecewiseSpace,n)
   k=div(n,length(d))
    r=n-length(d)*k

    [vcat([points(d.spaces[j],k+1) for j=1:r]...);
        vcat([points(d.spaces[j],k) for j=r+1:length(d)]...)]
end

function transform(S::PiecewiseSpace,vals::Vector,plan...)
    n=length(vals)
    K=length(S)
    k=div(n,K)
    PT=coefficient_type(S,eltype(vals))
    if k==0
        M=Array(Vector{PT},n)
        for j=1:n
            M[j]=transform(S[j],[vals[j]])
        end
    else
        r=n-K*k
        M=Array(Vector{PT},K)

        for j=1:r
            M[j]=transform(S[j],vals[(j-1)*(k+1)+1:j*(k+1)])
        end
        for j=r+1:length(S)
            M[j]=transform(S[j],vals[r*(k+1)+(j-r-1)*k+1:r*(k+1)+(j-r)*k])
        end
    end

    interlace(M,S)
end

itransform(S::PiecewiseSpace,cfs::Vector,plan...) =
    vcat([itransform(S.spaces[j],Fun(cfs,S)[j].coefficients) for j=1:length(S)]...)




itransform(S::SumSpace,cfs,plan...) = Fun(cfs,S)(points(S,length(cfs)))
itransform!(S::SumSpace,cfs,plan...) = (cfs[:]=Fun(cfs,S)(points(S,length(cfs))))



## SumSpace{ConstantSpace}
# this space is special

union_rule(P::PiecewiseSpace,C::ConstantSpace{AnyDomain})=PiecewiseSpace(map(sp->union(sp,C),P.spaces))



## choosedomainspace

function choosedomainspace{T}(A::InterlaceOperator{T,1},sp::TupleSpace)
    # this ensures correct dispatch for unino
    sps = Vector{Space}(
        filter(x->!isambiguous(x),map(choosedomainspace,A.ops,sp.spaces)))
    if isempty(sps)
        UnsetSpace()
    else
        union(sps...)
    end
end
