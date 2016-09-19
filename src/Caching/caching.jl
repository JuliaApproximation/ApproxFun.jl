#### Caching
function CachedOperator{T}(io::InterlaceOperator{T,1};padding::Bool=false)
    ds=domainspace(io)
    rs=rangespace(io)

    ind=find(op->isinf(size(op,1)),io.ops)
    if length(ind) ≠ 1  || !isbanded(io.ops[ind[1]])  # is almost banded
        return default_CachedOperator(io;padding=padding)
    end
    i=ind[1]
    bo=io.ops[i]
    lin,uin=bandwidths(bo)



    # calculate number of rows interlaced
    # each each row in the rangespace increases the lower bandwidth by 1
    nds=0
    md=0
    for k=1:length(io.ops)
        if k ≠ i
            d=dimension(rs[k])
            nds+=d
            md=max(md,d)
        end
    end


    isend=true
    for k=i+1:length(io.ops)
        if dimension(rs[k]) == md
            isend=false
        end
    end

    numoprows=isend?md-1:md
    n=nds+numoprows

    (l,u) = (max(lin+nds,n-1),max(0,uin+1-ind[1]))

    # add extra rows for QR
    if padding
        u+=l
    end


    ret=abzeros(T,n,n+u,l,u,nds)

    # populate the finite rows
    jr=1:n+u
    ioM=io[1:n,jr]


    bcrow=1
    oprow=0
    for k=1:n
        K,J=io.rangeinterlacer[k]

        if K ≠ i
            # fill the fill matrix
            ret.fill.V[:,bcrow] = Matrix(view(io.ops[K],J:J,jr))
            ret.fill.U[k,bcrow] = 1
            bcrow += 1
        else
            oprow+=1
        end


        for j=rowrange(ret.bands,k)
            ret[k,j] = ioM[k,j]
        end
    end


    CachedOperator(io,ret,(n,n+u),ds,rs,(l,∞))
end


function CachedOperator{T}(io::InterlaceOperator{T,2};padding::Bool=false)
    ds=domainspace(io)
    rs=rangespace(io)
    di=io.domaininterlacer
    ri=io.rangeinterlacer
    ddims=dimensions(di.iterator)
    rdims=dimensions(ri.iterator)

    # we are only almost banded if every operator is either finite
    # range or banded, and if the # of  ∞ spaces is the same
    # for the domain and range
    isab=all(op->isfinite(size(op,1)) || isbanded(op),io.ops) &&
        count(isinf,ddims) == count(isinf,rdims)

    if !isab
        return CachedOperator(io,Array(T,0,0))
    end


    # these are the bandwidths if we only had ∞-dimensional operators
    d∞=find(isinf,[ddims...])
    r∞=find(isinf,[rdims...])
    p=length(d∞)

    for k in d∞
        @assert blocklengths(ds[k]) == Repeated(true)
    end
    for k in r∞
        @assert blocklengths(rs[k]) == Repeated(true)
    end

    l∞,u∞ = 0,0
    for k=1:p,j=1:p
        l∞=min(l∞,p*bandinds(io.ops[r∞[k],d∞[j]],1)+j-k)
    end
    for k=1:p,j=1:p
        u∞=max(u∞,p*bandinds(io.ops[r∞[k],d∞[j]],2)+j-k)
    end

    # now we move everything by the finite rank
    ncols=mapreduce(d->isfinite(d)?d:0,+,ddims)
    nbcs=mapreduce(d->isfinite(d)?d:0,+,rdims)
    shft=ncols-nbcs
    l∞,u∞=l∞+shft,u∞+shft

    # iterate through finite rows to find worst case bandinds
    l,u=l∞,u∞
        for k=1+nbcs+p,j=1:ncols+p
            N,n=ri[k]
            M,m=di[j]
            l=min(l,bandinds(io.ops[N,M],1)+n-m-k+j)
            u=max(u,bandinds(io.ops[N,M],2)+n-m-k+j)
        end


    # add extra rows for QR
    if padding
        u-=l
    end

    n=1+nbcs+p
    ret=abzeros(T,n,n+u,-l,u,nbcs)

    # populate entries and fill functionals
    bcrow=1
    oprow=0
    jr=1:n+u
    for k=1:n
        K,J=io.rangeinterlacer[k]

        if isfinite(rdims[K] )
            # fill the fill matrix
            ret.fill.V[:,bcrow] = Matrix(view(io,k:k,jr))
            ret.fill.U[k,bcrow] = 1
            bcrow += 1
        else
            oprow+=1
        end


        for j=rowrange(ret.bands,k)
            ret[k,j] = io[k,j]
        end
    end

    CachedOperator(io,ret,(n,n+u),ds,rs,(l,∞))
end



# back substitution
trtrs!(::Type{Val{'U'}},co::CachedOperator,u::Array) =
                trtrs!(Val{'U'},resizedata!(co,size(u,1),size(u,1)).data,u)




## Ac_mul_B! for QROperatorQ

function Ac_mul_Bpars(A::QROperatorQ,B::Vector,tolerance,maxlength)
    T = promote_type(eltype(A),eltype(B))
    Ac_mul_Bpars(Operator{T}(A),Vector{T}(B),tolerance,maxlength)
end


include("matrix.jl")
include("ragged.jl")
include("banded.jl")
include("almostbanded.jl")
include("bandedblock.jl")
