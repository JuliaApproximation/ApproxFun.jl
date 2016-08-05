export ⊕,depiece,pieces,PiecewiseSpace


## Interlace interator
# this interlaces the entries treating everything equally,
# but will stop interlacing finite dimensional when it runs
# out.
# Dimensions are either Int or ∞

immutable InterlaceIterator{DMS<:Tuple}
    dimensions::DMS
end

InterlaceIterator(V::Vector) = InterlaceIterator(tuple(V...))

Base.start(it::InterlaceIterator) = (findfirst(d-> d≠0,it.dimensions),1)
function Base.next(it::InterlaceIterator,st)
    if st[1] == length(it.dimensions)
        k = 1
        j = st[2]+1
    else
        k = st[1]+1
        j = st[2]
    end
    nst=(k,j)
    if done(it,nst)
        (st,nst)
    else
        # call next if we are not done
        (st,j > it.dimensions[k]?next(it,nst)[2]:nst)
    end
end

# are all Ints, so finite dimensional
function Base.done{N}(it::InterlaceIterator{NTuple{N,Int}},st)
    for k=1:length(it.dimensions)
        if st[2] ≤ it.dimensions[k] && k ≤ st[1]
            return false
        end
    end
    return true
end

# Are not all Ints, so ∞ dimensional
Base.done(it::InterlaceIterator,st) = false
Base.eltype(it::InterlaceIterator) = Tuple{Int,Int}





## SumSpace encodes a space that can be decoupled as f(x) = a(x) + b(x) where a is in S and b is in V


abstract DirectSumSpace{SV,T,DD,d} <: Space{T,DD,d}


dimension(sp::DirectSumSpace) = mapreduce(dimension,+,sp.spaces)

InterlaceIterator(sp::DirectSumSpace) = InterlaceIterator(map(dimension,sp.spaces))
interlacer(sp::DirectSumSpace) = InterlaceIterator(sp)
interlacer(sp::Space) = InterlaceIterator(tuple(dimension(sp)))
cache(Q::InterlaceIterator) = CachedIterator(Q)


for TYP in (:SumSpace,:TupleSpace)
    @eval begin
        immutable $TYP{SV,T,DD,d} <: DirectSumSpace{SV,T,DD,d}
            spaces::SV
            $TYP(dom::Domain) = new(tuple(map(typ->typ(dom),SV.parameters)...))
            $TYP(sp::Tuple) = new(sp)
        end

        $TYP(sp::Tuple) = $TYP{typeof(sp),mapreduce(basistype,promote_type,sp),
                             typeof(domain(first(sp))),ndims(first(sp))}(sp)
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
                               typeof(UnionDomain(map(domain,sp))),ndims(first(sp))}(sp)
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


for OP in (:(Base.getindex),:(Base.length),:(Base.next),:(Base.start),:(Base.done),:(Base.endof))
    @eval $OP(S::DirectSumSpace,k...)=$OP(S.spaces,k...)
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

Base.ones{T<:Number,SS,V}(::Type{T},S::PiecewiseSpace{SS,V}) = depiece(map(ones,S.spaces))
Base.ones(S::PiecewiseSpace)=ones(Float64,S)


identity_fun(S::PiecewiseSpace) = depiece(map(identity_fun,S.spaces))

# vec

function Base.getindex{DSS<:DirectSumSpace}(f::Fun{DSS},k::Integer)
    sp=space(f).spaces
    it=InterlaceIterator(map(dimension,sp))
    N=length(f.coefficients)
    d=dimension(sp[k])

    # preallocate: we know we have at most N coefficients
    ret=Array(eltype(f),N)
    j=1  # current coefficient
    p=0  # current length
    for (n,m) in it
        if j > N || m > d
            break
        end
        if n==k
            ret[m]=f.coefficients[j]
            p+=1
        end
        j+=1
    end
    resize!(ret,p)  # through out extra coefficients
    Fun(ret,sp[k])
end


# interlace coefficients according to iterator
function interlace{T,V<:AbstractVector}(::Type{T},v::AbstractVector{V},it::InterlaceIterator)
    ret=Vector{T}()
    N=mapreduce(length,max,v)
    for (n,m) in it
        if m > N
            break
        end

        if m ≤ length(v[n])
            push!(ret,v[n][m])
        else
            push!(ret,zero(T))
        end
    end
    ret
end

interlace{V<:AbstractVector}(v::AbstractVector{V},it::InterlaceIterator) =
    interlace(mapreduce(eltype,promote_type,v),v,it)

interlace{F<:Fun}(v::AbstractVector{F},sp::DirectSumSpace) =
    interlace(map(coefficients,v),InterlaceIterator(sp))

function interlace{F<:Fun}(v::Tuple{Vararg{F}},sp::DirectSumSpace)
    V=Array(Vector{mapreduce(eltype,promote_type,v)},length(v))
    for k=1:length(v)
        V[k]=coefficients(v[k])
    end
    interlace(V,InterlaceIterator(sp))
end


Base.vec(S::DirectSumSpace) = S.spaces
Base.vec{S<:DirectSumSpace}(f::Fun{S}) = Fun[f[j] for j=1:length(space(f).spaces)]

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

        $Dep(v::Vector{Any})=depiece([v...])
    end
end



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
        ret=Array(PT,n)
        for j=1:n
            ret[j]=transform(S[j],[vals[j]])[1]
        end

        ret
    else
        r=n-K*k
        M=Array(PT,k+1,K)

        for j=1:r
            M[:,j]=transform(S[j],vals[(j-1)*(k+1)+1:j*(k+1)])
        end
        for j=r+1:length(S)
            M[1:k,j]=transform(S[j],vals[r*(k+1)+(j-r-1)*k+1:r*(k+1)+(j-r)*k])
            M[k+1,j]=zero(PT)
        end

    vec(M.')
    end
end

itransform(S::PiecewiseSpace,cfs::Vector,plan...)=vcat([itransform(S.spaces[j],cfs[j:length(S):end]) for j=1:length(S)]...)




itransform(S::SumSpace,cfs,plan...)=Fun(cfs,S)(points(S,length(cfs)))
itransform!(S::SumSpace,cfs,plan...)=(cfs[:]=Fun(cfs,S)(points(S,length(cfs))))



## SumSpace{ConstantSpace}
# this space is special

union_rule(P::PiecewiseSpace,C::ConstantSpace{AnyDomain})=PiecewiseSpace(map(sp->union(sp,C),P.spaces))
