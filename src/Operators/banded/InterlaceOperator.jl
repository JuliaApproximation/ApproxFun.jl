



##interlace block operators

iszerooperator(::ZeroOperator) = true
iszerooperator(A::ConstantOperator) = A.c==0.
iszerooperator(A) = false
function isboundaryrow(A,k)
    for j=1:size(A,2)
        if isafunctional(A[k,j])
            return true
        end
    end

    return false
end



domainscompatible{T<:Operator}(A::Matrix{T})=domainscompatible(map(domainspace,A))

function spacescompatible{T<:Operator}(A::Matrix{T})
    for k=1:size(A,1)
        if !spacescompatible(map(rangespace,vec(A[k,:])))
            return false
        end
    end
    for k=1:size(A,2)
        if !spacescompatible(map(domainspace,A[:,k]))
            return false
        end
    end
    true
end

spacescompatible{T<:Operator}(A::Vector{T}) = spacescompatible(map(domainspace,A))

function domainspace{T<:Operator}(A::Matrix{T})
    @assert spacescompatible(A)

    spl=map(domainspace,vec(A[1,:]))
    if spacescompatible(spl)
        ArraySpace(first(spl),length(spl))
    elseif domainscompatible(spl)
        TupleSpace(spl)
    else
        PiecewiseSpace(spl)
    end
end

function rangespace{T<:Operator}(A::Vector{T})
    @assert spacescompatible(A)

    spl=map(rangespace,A)
    if spacescompatible(spl)
        ArraySpace(first(spl),length(spl))
    elseif domainscompatible(spl)
        TupleSpace(spl)
    else
        PiecewiseSpace(spl)
    end
end

function promotespaces{T<:Operator}(A::Matrix{T})
    A=copy(A)#TODO: promote might have different Array type
    for j=1:size(A,2)
        A[:,j]=promotedomainspace(A[:,j])
    end
    for k=1:size(A,1)
        A[k,:]=promoterangespace(vec(A[k,:]))
    end

    # do a second loop as spaces might have been inferred
    # during range space
    for j=1:size(A,2)
        A[:,j]=promotedomainspace(A[:,j])
    end
    A
end


## Interlace operator

immutable InterlaceOperator{T,DS,RS,DI,RI,BI} <: Operator{T}
    ops::Matrix{Operator{T}}
    domainspace::DS
    rangespace::RS
    domaininterlacer::DI
    rangeinterlacer::RI
    bandinds::BI
end
function InterlaceOperator{T}(ops::Matrix{Operator{T}},ds::Space,rs::Space)
    # calculate bandinds
    p=size(ops,1)
    if size(ops,2) == p && all(isbanded,ops)
        l,u = 0,0
        for k=1:p,j=1:p
            l=min(l,p*bandinds(ops[k,j],1)+j-k)
        end
        for k=1:p,j=1:p
            u=max(u,p*bandinds(ops[k,j],2)+j-k)
        end
    else
        l,u = (-∞,∞)  # not banded
    end


    InterlaceOperator(ops,ds,rs,
                        CachedIterator(interlacer(ds)),
                        CachedIterator(interlacer(rs)),
                        (l,u))
end

function InterlaceOperator{T}(opsin::Matrix{Operator{T}})
    ops=promotespaces(opsin)
    # TODO: make consistent
    # if its a row vector, we assume scalar
    if size(ops,1) == 1
        InterlaceOperator(ops,domainspace(ops),rangespace(ops[1]))
    else
        InterlaceOperator(ops,domainspace(ops),rangespace(ops[:,1]))
    end
end

function InterlaceOperator{T}(opsin::Vector{Operator{T}})
    ops=promotedomainspace(opsin)
    InterlaceOperator(ops'',domainspace(first(ops)),rangespace(ops))
end

InterlaceOperator(ops::Matrix) =
    InterlaceOperator(Matrix{Operator{mapreduce(eltype,promote_type,ops)}}(ops))

function Base.convert{T}(::Type{Operator{T}},S::InterlaceOperator)
    if T == eltype(S)
        S
    else
        ops=Array(Operator{T},size(S.ops)...)
        for j=1:size(S.ops,2),k=1:size(S.ops,1)
            ops[k,j]=S.ops[k,j]
        end
        InterlaceOperator(ops,domainspace(S),rangespace(S),
                            S.domaininterlacer,S.rangeinterlacer,S.bandinds)
    end
end



#TODO: More efficient to save bandinds
bandinds(M::InterlaceOperator) = M.bandinds

function getindex{T}(op::InterlaceOperator{T},k::Integer,j::Integer)
    M,J = op.domaininterlacer[j]
    N,K = op.rangeinterlacer[k]
    op.ops[N,M][K,J]::T
end

function getindex{T}(op::InterlaceOperator{T},k::Integer)
    if size(op,1) == 1
        op[1,k]
    elseif size(op,2) == 1
        op[k,1]
    else
        error("Only implemented for row/column operators.")
    end
end


domainspace(IO::InterlaceOperator) = IO.domainspace
rangespace(IO::InterlaceOperator) = IO.rangespace

#tests whether an operator can be made into a column
iscolop(op) = isconstop(op)
iscolop(::Multiplication) = true


function interlace{T<:Operator}(A::Matrix{T})
    if all(isbanded,A)
        # Hack to use BlockOperator when appropriate
        if size(A,1)==1 && all(iscolop,A[1,1:end-1])
            return BlockOperator(A)
        else
            return InterlaceOperator(A)
        end
    end


    m,n=size(A)
    TT=mapreduce(eltype,promote_type,A)
    # Use BlockOperator whenever the first columns are all constants
    if all(isconstop,A[1:end-1,1:end-1]) &&
            all(iscolop,A[end,1:end-1]) &&
            all(a->isafunctional(a),A[1:end-1,end]) && isbanded(A[end,end]) &&
            isinf(size(A[end],1)) && isinf(size(A[end],2))
        return [Operator{TT}[BlockFunctional(map(Number,vec(A[k,1:end-1])),A[k,end]) for k=1:m-1]...;
                    BlockOperator(A[end:end,:])]
    end


    A=promotespaces(A)

    br=0#num boundary rows
    for k=1:m
        if isboundaryrow(A,k)
            br+=1
        end
    end

    for k=1:br
        @assert isboundaryrow(A,k)
    end

    S=Array(Operator{TT},br<m?br+1:br)

    for k=1:br
        S[k] = InterlaceOperator(A[k,:])
    end

    if br < m
        Am=A[br+1:m,:]
        S[br+1] = interlace(Am)
    end

    if(size(S,1) ==1)
        S[1]
    else
        S
    end
end



immutable DiagonalInterlaceOperator{OPS,DS,RS,T<:Number} <: Operator{T}
    ops::OPS
    domainspace::DS
    rangespace::RS
end

function DiagonalInterlaceOperator(v::Tuple,ds::Space,rs::Space)
    T=mapreduce(eltype,promote_type,v)
    w=map(Operator{T},v)
    DiagonalInterlaceOperator{typeof(w),typeof(ds),typeof(rs),T}(w,ds,rs)
end
DiagonalInterlaceOperator{ST<:Space}(v::Tuple,::Type{ST})=DiagonalInterlaceOperator(v,ST(map(domainspace,v)),ST(map(rangespace,v)))
DiagonalInterlaceOperator(v::Vector,k...)=DiagonalInterlaceOperator(tuple(v...),k...)


Base.convert{T}(::Type{Operator{T}},op::DiagonalInterlaceOperator)=
        DiagonalInterlaceOperator(map(Operator{T},op.ops),op.domainspace,op.rangespace)


function bandinds(S::DiagonalInterlaceOperator)
    binds=map(bandinds,S.ops)
    bra=mapreduce(first,min,binds)
    brb=mapreduce(last,max,binds)
    n=length(S.ops)
    n*bra,n*brb
end


function getindex(D::DiagonalInterlaceOperator,k::Integer,j::Integer)
    n=length(D.ops)
    mk = n+mod(k,-n)
    T=eltype(D)
    if mk == n+mod(j,-n)  # same block
        k=(k-1)÷n+1  # map k and j to block coordinates
        j=(j-1)÷n+1
        D.ops[mk][k,j]::T
    else
        zero(T)
    end
end

domainspace(D::DiagonalInterlaceOperator)=D.domainspace
rangespace(D::DiagonalInterlaceOperator)=D.rangespace
