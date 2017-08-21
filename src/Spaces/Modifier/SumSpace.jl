export ⊕,components,PiecewiseSpace



## BlockInterlacer
# interlaces coefficients by blocks
# this has the property that all the coefficients of a block of a subspace
# are grouped together, starting with the first bloc
#
# TODO: cache sums


struct BlockInterlacer{DMS<:Tuple}
    blocks::DMS
end


const TrivialInterlacer{d} = BlockInterlacer{NTuple{d,Repeated{Bool}}}

BlockInterlacer(v::AbstractVector) = BlockInterlacer(tuple(v...))

Base.eltype(it::BlockInterlacer) = Tuple{Int,Int}

dimensions(b::BlockInterlacer) = map(sum,b.blocks)
dimension(b::BlockInterlacer,k) = sum(b.blocks[k])
Base.length(b::BlockInterlacer) = mapreduce(sum,+,b.blocks)

# the state is always (whichblock,curblock,cursubblock,curcoefficients)
Base.start(it::BlockInterlacer) = (1,1,map(start,it.blocks),ntuple(zero,length(it.blocks)))

function Base.next(it::BlockInterlacer,st)
    N,k,blkst,lngs = st

    if N>length(it.blocks)
        # increment to next block
        blkst = map((blit,blst)->done(blit,blst) ? blst : next(blit,blst)[2],it.blocks,blkst)
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
        if st[end][k] < sum(it.blocks[k])
            return false
        end
    end
    return true
end






## SumSpace encodes a space that can be decoupled as f(x) = a(x) + b(x) where a is in S and b is in V


abstract type DirectSumSpace{SV,D,R} <: Space{D,R} end


dimension(sp::DirectSumSpace) = mapreduce(dimension,+,sp.spaces)

components(s::Space) = (s,)
components(sp::DirectSumSpace) = sp.spaces
component(s::Space,k...) = components(s)[k...]

ncomponents(s::Space) = length(components(s))
ncomponents(f::Fun) = ncomponents(space(f))

BlockInterlacer(sp::DirectSumSpace) = BlockInterlacer(map(blocklengths,sp.spaces))
interlacer(sp::DirectSumSpace) = BlockInterlacer(sp)
interlacer(sp::Space) = BlockInterlacer(tuple(blocklengths(sp)))
cache(Q::BlockInterlacer) = CachedIterator(Q)

function blocklengths(sp::DirectSumSpace)
    bl=map(blocklengths,components(sp))
    N=mapreduce(length,max,bl)
    mapreduce(b->pad(b,N),+,bl)
end
block(sp::DirectSumSpace,k::Int)::Block = findfirst(x->x≥k,cumsum(blocklengths(sp)))



isambiguous(sp::DirectSumSpace) = any(isambiguous,components(sp))


struct SumSpace{SV,D,R} <: DirectSumSpace{SV,D,R}
    spaces::SV
    SumSpace{SV,D,R}(dom::Domain) where {SV,D,R} =
        new(tuple(map(typ->typ(dom),SV.parameters)...))
    SumSpace{SV,D,R}(sp::Tuple) where {SV,D,R} = new(sp)
end

SumSpace(sp::Tuple) = SumSpace{typeof(sp),domaintype(first(sp)),
                                mapreduce(rangetype,promote_type,sp)}(sp)


struct PiecewiseSpace{SV,D<:UnionDomain,R} <: DirectSumSpace{SV,D,R}
    spaces::SV
    PiecewiseSpace{SV,D,R}(dom::AnyDomain) where {SV,D,R} =
        new{SV,D,R}(tuple(map(typ->typ(dom),SV.parameters)...))
    PiecewiseSpace{SV,D,R}(dom::UnionDomain) where {SV,D,R} =
        new{SV,D,R}(tuple(map((typ,dom)->typ(dom),SV.parameters,dom.domains)...))
    PiecewiseSpace{SV,D,R}(sp::Tuple) where {SV,D,R} =
        new{SV,D,R}(sp)
end

function PiecewiseSpace(spin::Tuple)
    sp=tuple(union(spin)...)  # remove duplicates

    PiecewiseSpace{typeof(sp),typeof(UnionDomain(map(domain,sp))),
                   mapreduce(rangetype,promote_type,sp)}(sp)
end



for TYP in (:SumSpace,:PiecewiseSpace)
    @eval begin
        $TYP(A::$TYP,B::$TYP) = $TYP(tuple(A.spaces...,B.spaces...))

        $TYP(A::Space,B::$TYP) = $TYP(tuple(A,B.spaces...))
        $TYP(A::$TYP,B::Space) = $TYP(tuple(A.spaces...,B))
        $TYP(A::Space...) = $TYP(A)
        $TYP(sp::AbstractArray) = $TYP(tuple(sp...))

        canonicalspace(A::$TYP) = $TYP(sort([A.spaces...]))
    end
end



setdomain(A::SumSpace,d::Domain) = SumSpace(map(sp->setdomain(sp,d),A.spaces))


setdomain(A::PiecewiseSpace,d::UnionDomain) =
    PiecewiseSpace(map((sp,dd)->setdomain(sp,dd),A.spaces,d.domains))



function spacescompatible(A::S,B::S) where S<:DirectSumSpace
    if ncomponents(A) != ncomponents(B)
        false
    else
        ret=true
        for k=1:ncomponents(A)
            if !spacescompatible(component(A,k),component(B,k))
                return false
            end
        end
        true
    end
end

function Base.promote_rule(::Type{Fun{SumSpace{SV,D,R},V,VV}},::Type{T}) where {SV,D,R,V,T<:Number,VV}
    for k=1:length(SV.parameters)
        pt=promote_type(VFun{SV.parameters[k],V},T)
        if pt != Fun
            return VFun{SumSpace{Tuple{SV.parameters[1:k-1]...,pt.parameters[1],SV.parameters[k+1:end]...},
                       D,R},promote_type(V,T)}
        end
    end
    Fun
end

Base.promote_rule(::Type{Fun{SumSpace{SV,D,R}}},::Type{T}) where {SV,D,R,T<:Number} =
    promote_rule(VFun{SumSpace{SV,D,R},Float64},T)

function Base.promote_rule(::Type{Fun{PiecewiseSpace{SV,D,R},V,VV}},::Type{T}) where {SV,D,R,V,VV,T<:Number}
    # if any doesn't support promoting, just leave unpromoted

    newfsp=map(s->promote_type(VFun{s,V},T),SV.parameters)
    if any(s->s==Fun,newfsp)
        Fun
    else
        newsp=map(s->s.parameters[1],newfsp)
        VFun{PiecewiseSpace{Tuple{newsp...},D,R},promote_type(V,T)}
    end
end

Base.promote_rule(::Type{Fun{PiecewiseSpace{SV,D,R}}},::Type{T}) where {SV,D,R,T<:Number} =
    promote_rule(VFun{PiecewiseSpace{SV,D,R},Float64},T)



# domain
domain(A::SumSpace) = domain(A.spaces[end])      # TODO: this assumes all spaces have the same domain
                                                     #        we use end to avoid ConstantSpace


Space(d::UnionDomain) = PiecewiseSpace(map(Space,d.domains))
domain(S::PiecewiseSpace) = UnionDomain(map(domain,S.spaces))



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


## components

# We use a view when it's avaiable to avoid allocation
component_coefficients(it::TrivialInterlacer{d},cfs,k) where {d} = (@view cfs[k:d:end])

function component_coefficients(it,cfs,k)
    N=length(cfs)
    d=dimension(it,k)

    # preallocate: we know we have at most N coefficients
    ret=Array{eltype(cfs)}(N)
    j=1  # current coefficient
    p=0  # current length
    for (n,m) in it
        if j > N
            break
        end
        if n==k
            ret[m] = cfs[j]
            p+=1
            if m ≥ d
                # if we've reached the dimension, we are done
                break
            end
        end
        j+=1
    end
    resize!(ret,p)  # throw out extra coefficients
end

component_coefficients(sp::Space,cfs,k) = component_coefficients(interlacer(sp),cfs,k)

component(f::Fun{<:DirectSumSpace},k::Integer) =
    Fun(component(space(f),k),component_coefficients(space(f),f.coefficients,k))


## evaluate


for OP in (:(Base.last),:(Base.first))
    @eval begin
        $OP(f::Fun{SS}) where {SS<:SumSpace} = mapreduce($OP,+,components(f))
        $OP(f::Fun{SS}) where {SS<:PiecewiseSpace} = $OP($OP(components(f)))
    end
end

# this is a type-stable version
# TODO: replace with generated function
function evaluate(f::AbstractVector,S::SumSpace{Tuple{A,B}},x) where {A,B}
    it = interlacer(S)
    a,b = S.spaces
    evaluate(component_coefficients(it,f,1),a,x) +
        evaluate(component_coefficients(it,f,2),b,x)
end

function evaluate(f::AbstractVector,S::SumSpace,x)
    ret = zero(rangetype(S))
    it = interlacer(S)
    for k=1:ncomponents(S)
        ret += evaluate(component_coefficients(it,f,k),component(S,k),x)
    end
    ret
end


function evaluate(f::AbstractVector,S::PiecewiseSpace,x)
    d=domain(S)
    g=Fun(S,f)

#    ret=zero(promote_type(eltype(f),eltype(S)))
    for k=1:ncomponents(d)
        if in(x,component(d,k))
            return component(g,k)(x)
        end
    end
    return 0*first(g)
end


## calculus
for TYP in (:SumSpace,:PiecewiseSpace)
    for OP in (:differentiate,:integrate)
        @eval function $OP(f::Fun{D,T}) where {D<:$TYP,T}
            fs = map($OP,components(f))
            sp = $TYP(map(space,fs))
            Fun(sp,interlace(fs,sp))
        end
    end
    for OP in (:(Base.real),:(Base.imag),:(Base.conj))
        @eval begin
            $OP(f::Fun{$TYP{SV,DD,RR}}) where {SV,DD,RR<:Real} = Fun(f.space,$OP(f.coefficients))
            function $OP(f::Fun{$TYP{SV,DD,RR}}) where {SV,DD,RR}
                fs=map($OP,components(f))
                sp=$TYP(map(space,fs))
                Fun(sp,interlace(fs,sp))
            end
        end
    end
end



@eval function Base.cumsum(f::Fun{D,T}) where {D<:SumSpace,T}
    fs=map(cumsum,components(f))
    sp=SumSpace(map(space,fs))
    Fun(sp,interlace(fs,sp))
end


for TYP in (:SumSpace,:PiecewiseSpace), OP in (:(Base.sum),:linesum)
    @eval $OP(f::Fun{V}) where {V<:$TYP} = mapreduce($OP,+,components(f))
end

function Base.cumsum(f::Fun{V}) where V<:PiecewiseSpace
    vf=components(f)
    r=zero(eltype(f))
    for k=1:length(vf)
        vf[k]=cumsum(vf[k]) + r
        r=last(vf[k])
    end
    Fun(vf,PiecewiseSpace)
end

Base.cumsum(f::Fun{V},d::Domain) where {V<:PiecewiseSpace} =
    mapreduce(g->cumsum(g,d),+,components(f))



bilinearform(f::Fun{S},g::Fun{V}) where {S<:PiecewiseSpace,V<:PiecewiseSpace} =
    sum(map(bilinearform,components(f),components(g)))
linebilinearform(f::Fun{S},g::Fun{V}) where {S<:PiecewiseSpace,V<:PiecewiseSpace} =
    sum(map(linebilinearform,components(f),components(g)))

# assume first domain has 1 as a basis element



function Base.ones(::Type{T},S::SumSpace) where T<:Number
    @assert ncomponents(S) == 2
    if isconvertible(ConstantSpace(),component(S,1))
        ones(T,component(S,1)) ⊕ zeros(T,component(S,2))
    else
        zeros(T,component(S,1)) ⊕ ones(T,component(S,2))
    end
end

Base.ones(S::SumSpace) = ones(Float64,S)

Base.ones(::Type{T},S::PiecewiseSpace{SS,V}) where {T<:Number,SS,V} =
    Fun(map(ones,components(S)),PiecewiseSpace)
Base.ones(S::PiecewiseSpace) = ones(Float64,S)


identity_fun(S::PiecewiseSpace) = Fun(map(identity_fun,S.spaces),PiecewiseSpace)


# interlace coefficients according to iterator
function interlace(::Type{T},v::AbstractVector{V},it::BlockInterlacer) where {T,V<:AbstractVector}
    ret=Array{T}(0)
    N=mapreduce(length,max,v)
    cnts = Vector(map(length,v))  # convert to Vector to ensure mutable

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

interlace(v::AbstractVector{V},it::BlockInterlacer) where {V<:AbstractVector} =
    interlace(mapreduce(eltype,promote_type,v),v,it)

interlace(v::AbstractVector{V},sp::DirectSumSpace) where {V<:AbstractVector} =
    interlace(v,interlacer(sp))

interlace(v::AbstractVector{F},sp::DirectSumSpace) where {F<:Fun} =
    interlace(map(coefficients,v),sp)


function interlace(v::Union{Tuple,Vector{Any}},sp::DirectSumSpace)
    if all(vk->isa(vk,Fun),v)
        V=Array{Vector{mapreduce(eltype,promote_type,v)}}(length(v))
        for k=1:length(v)
            V[k]=coefficients(v[k])
        end
    else
        V=Array{Vector{mapreduce(eltype,promote_type,v)}}(length(v))
        for k=1:length(v)
            V[k]=v[k]
        end
    end
    interlace(V,sp)
end

components(f::Fun{S}) where {S<:DirectSumSpace} = Fun[component(f,j) for j=1:ncomponents(f)]

function Fun(v::AbstractVector{F},::Type{PiecewiseSpace}) where F<:Fun
    sp = PiecewiseSpace(map(space,v))
    Fun(sp,interlace(v,sp))
end
function Fun(v::Tuple,::Type{PiecewiseSpace})
    sp=PiecewiseSpace(map(space,v))
    Fun(sp,interlace(v,sp))
end

Fun(v::AbstractVector{Any},::Type{PiecewiseSpace}) = Fun(tuple(v...),PiecewiseSpace)

## transforms


function points(d::PiecewiseSpace,n)
   k=div(n,ncomponents(d))
    r=n-ncomponents(d)*k

    [vcat([points(d.spaces[j],k+1) for j=1:r]...);
        vcat([points(d.spaces[j],k) for j=r+1:ncomponents(d)]...)]
end

plan_transform(sp::PiecewiseSpace,vals::AbstractVector) =
    TransformPlan{eltype(vals),typeof(sp),false,Void}(sp,nothing)

plan_itransform(sp::PiecewiseSpace,vals::AbstractVector) =
    ITransformPlan{eltype(vals),typeof(sp),false,Void}(sp,nothing)



function *(P::TransformPlan{T,PS,false},vals::AbstractVector{T}) where {PS<:PiecewiseSpace,T}
    S=components(P.space)
    n=length(vals)
    K=length(S)
    k=div(n,K)
    PT=promote_type(prectype(P.space),eltype(vals))
    if k==0
        M=Array{Vector{PT}}(n)
        for j=1:n
            M[j]=transform(S[j],[vals[j]])
        end
    else
        r=n-K*k
        M=Array{Vector{PT}}(K)

        for j=1:r
            M[j]=transform(S[j],vals[(j-1)*(k+1)+1:j*(k+1)])
        end
        for j=r+1:length(S)
            M[j]=transform(S[j],vals[r*(k+1)+(j-r-1)*k+1:r*(k+1)+(j-r)*k])
        end
    end

    interlace(M,P.space)
end

*(P::ITransformPlan{T,PS,false},cfs::AbstractVector{T}) where {T,PS<:PiecewiseSpace} =
    vcat([itransform(P.space.spaces[j],component(Fun(P.space,cfs),j).coefficients) for j=1:ncomponents(P.space)]...)



itransform(S::SumSpace,cfs) = Fun(S,cfs).(points(S,length(cfs)))
function itransform!(S::SumSpace,cfs)
    f    = Fun(S,cfs)
    cfs .= f.(points(S,length(cfs)))
end



## SumSpace{ConstantSpace}
# this space is special

union_rule(P::PiecewiseSpace,C::ConstantSpace{AnyDomain}) =
    PiecewiseSpace(map(sp->union(sp,C),P.spaces))



## Multiplication
function *(f::Fun{<:PiecewiseSpace,T},g::Fun{<:PiecewiseSpace,N}) where {T,N}
    domain(f) ≠ domain(g) && return default_mult(f,g)

    Fun(map(*,components(f),components(g)),PiecewiseSpace)
end



## Multivariate

ncomponents(sp::TensorSpace) = mapreduce(s -> ncomponents(s), *, factors(sp))

component(sp::TensorSpace{Tuple{S1,S2}},k::Integer) where {S1<:DirectSumSpace,S2<:DirectSumSpace} =
    error("Not defined. Used component(sp,k,j).")

component(sp::TensorSpace{Tuple{S1,S2}},k::Integer) where {S1<:DirectSumSpace,S2} =
    component(factor(sp,1),k) ⊗ factor(sp,2)

component(sp::TensorSpace{Tuple{S1,S2}},k::Integer) where {S1,S2<:DirectSumSpace} =
    factor(sp,1) ⊗ component(factor(sp,2),k)


component(sp::TensorSpace{Tuple{S1,S2}},k::Integer,j::Integer) where {S1<:DirectSumSpace,S2<:DirectSumSpace} =
    component(factor(sp,1),k) ⊗ component(factor(sp,2),j)
