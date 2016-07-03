

export StrideOperator,StrideFunctional



#S[rowstride*k + rowindex,colstride*j + colindex] == op[k,j]
#S[k,j] == op[(k-rowindex)/rowstride,(j-colindex)/colstride]

immutable StrideOperator{T<:Number,B<:Operator} <: BandedOperator{T}
    op::B
    rowindex::Int
    colindex::Int
    rowstride::Int
    colstride::Int

    function StrideOperator(o,r,c,rs,cs)
        @assert rs == cs
        @assert rs != 0

        new(o,r,c,rs,cs)
    end
end

StrideOperator{T<:Number}(B::Operator{T},r,c,rs,cs)=StrideOperator{T,typeof(B)}(B,r,c,rs,cs)
StrideOperator{T<:Number}(B::Operator{T},r,c,rs)=StrideOperator{T,typeof(B)}(B,r,c,rs,rs)
StrideOperator{T<:Number}(B::Operator{T},r,c)=StrideOperator{T,typeof(B)}(B,r,c,1,1)

function bandinds(S::StrideOperator)
    br=bandinds(S.op)

    st = abs(S.colstride)

    min(st*br[1]-S.rowindex+S.colindex,0),max(st*br[end]-S.rowindex+S.colindex,0)
end

Base.stride(S::StrideOperator)=mod(S.rowindex,S.rowstride)==mod(S.colindex,S.colstride)==0?S.rowstride:1



# First index above
firstrw(rs,ri,k::Integer)=fld(k-ri+rs-1,rs)
firstrw(S,k::Integer)=firstrw(S.rowstride,S.rowindex,k)

#Last index below
lastrw(rs,ri,k::Integer)=fld(k-ri,rs)


divrowrange(rs,ri,r)=max(1,firstrw(rs,ri,r[1])):lastrw(rs,ri,r[end])

for op in (:firstrw,:lastrw,:divrowrange)
    @eval $op(S,k...)=$op(S.rowstride,S.rowindex,k...)
end


#S[rowstride*k + rowindex,colstride*j + colindex] == op[k,j]
#S[k,j] == A[k,j-k]
#A[rowstride*k + rowindex,colstride*j + colindex - k] == op[k,j]

stride_getindex(op,ri,ci,rs,cs,k::Integer,j::Integer) =
    op[firstrw(rs,ri,k),firstrw(cs,ci,j)]



stride_getindex(S::StrideOperator,k::Integer,j::Integer) =
    stride_addentries!(S.op,S.rowindex,S.colindex,S.rowstride,S.colstride,k,j)


getindex(S::StrideOperator,k::Integer,j::Integer) = stride_addentries!(S,k,j)
domain(S::StrideOperator) = Any ##TODO: tensor product




## StrideFunctional


type StrideFunctional{T<:Number,B<:Operator} <: Operator{T}
    op::B
    rowindex::Int
    stride::Int
end

@functional StrideFunctional

StrideFunctional{T<:Number}(B::Operator{T},r,rs)=StrideFunctional{T,typeof(B)}(B,r,rs)


Base.getindex{T<:Number}(op::StrideFunctional{T},kr::Range)=T[((k-op.rowindex)≥1 &&(k-op.rowindex)%op.stride==0)?op.op[fld(k-op.rowindex,op.stride)]:zero(T) for k=kr]

@eval Base.convert{T}(::Type{Operator{T}},S::StrideFunctional)=StrideFunctional(convert(Operator{T},S.op),S.rowindex,S.stride)


##interlace block operators

iszerooperator(::ZeroOperator)=true
iszerooperator(::ZeroFunctional)=true
iszerooperator(A::ConstantOperator)=A.c==0.
iszerooperator(A)=false
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

spacescompatible{T<:Operator}(A::Vector{T})=spacescompatible(map(domainspace,A))

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

immutable InterlaceOperator{T} <: BandedOperator{T}
    ops::Matrix{BandedOperator{T}}
    function InterlaceOperator(os)
        @assert size(os,1)==size(os,2)
        new(promotespaces(os))
    end
end
InterlaceOperator{T}(ops::Matrix{BandedOperator{T}})=InterlaceOperator{T}(ops)
InterlaceOperator{B<:BandedOperator}(ops::Matrix{B})=InterlaceOperator(convert(Matrix{BandedOperator{mapreduce(eltype,promote_type,ops)}},ops))

function Base.convert{BT<:Operator}(::Type{BT},S::InterlaceOperator)
    ops=Array(BandedOperator{eltype(BT)},size(S.ops)...)
    for j=1:size(S.ops,2),k=1:size(S.ops,1)
        ops[k,j]=S.ops[k,j]
    end
    InterlaceOperator(ops)
end



#TODO: More efficient to save bandinds
function bandinds(M::InterlaceOperator,N::Integer)
    if N==1
        p=size(M.ops,1)
        ret=0
        for k=1:p,j=1:p
            ret=min(ret,p*bandinds(M.ops[k,j],1)+j-k)
        end
        ret
    else #N==2
        p=size(M.ops,1)
        ret=0
        for k=1:p,j=1:p
            ret=max(ret,p*bandinds(M.ops[k,j],2)+j-k)
        end
        ret
    end
end
bandinds(M::InterlaceOperator)=bandinds(M,1),bandinds(M,2)

function getindex(M::InterlaceOperator,k::Integer,j::Integer)
    n,m=size(M.ops)
    mk = n+mod(k,-n)
    mj = m+mod(j,-m)
    T=eltype(M)

    k=(k-1)÷n+1  # map k and j to block coordinates
    j=(j-1)÷m+1
    M.ops[mk,mj][k,j]::T
end

domainspace(IO::InterlaceOperator)=domainspace(IO.ops)
rangespace(IO::InterlaceOperator)=rangespace(IO.ops[:,1])

#tests whether an operator can be made into a column
iscolop(op)=isconstop(op)
iscolop(::Multiplication)=true

function interlace{T<:BandedOperator}(A::Matrix{T})
    # Hack to use BlockOperator when appropriate
    if size(A,1)==1 && all(iscolop,A[1,1:end-1])
        BlockOperator(A)
    else
        InterlaceOperator(A)
    end
end


# If the matrix is Operators, we assume it may contain
# functionals as well
function interlace{T<:Operator}(A::Matrix{T})
    m,n=size(A)
    TT=mapreduce(eltype,promote_type,A)
    # Use BlockOperator whenever the first columns are all constants
    if all(isconstop,A[1:end-1,1:end-1]) &&
            all(iscolop,A[end,1:end-1]) &&
            all(a->isafunctional(a),A[1:end-1,end]) && isa(A[end,end],BandedOperator)
        return [Operator{TT}[BlockFunctional(map(Number,vec(A[k,1:end-1])),A[k,end]) for k=1:m-1]...;
                    BlockOperator(A[end:end,:])]
    end


    A=promotespaces(A)

    dsp=domainspace(A)

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

    for k=1:br, j=1:n
        if !iszerooperator(A[k,j])
            op = StrideFunctional(A[k,j],j-n,n)

            if !isdefined(S,k)
                S[k] = op
            else
                S[k] = S[k] + op
            end
        end
    end

    for k=1:br
        S[k]=promotedomainspace(S[k],dsp)
    end

    if br < m
        Am=A[br+1:m,:]
        S[br+1]=interlace(convert(Matrix{BandedOperator{TT}},Am))
    end

    if(size(S,1) ==1)
        S[1]
    else
        S
    end
end



immutable DiagonalInterlaceOperator{OPS,DS,RS,T<:Number} <: BandedOperator{T}
    ops::OPS
    domainspace::DS
    rangespace::RS
end

DiagonalInterlaceOperator(v::Tuple,ds::Space,rs::Space)=DiagonalInterlaceOperator{typeof(v),typeof(ds),
                                                                                  typeof(rs),mapreduce(eltype,promote_type,v)}(v,ds,rs)
DiagonalInterlaceOperator{ST<:Space}(v::Tuple,::Type{ST})=DiagonalInterlaceOperator(v,ST(map(domainspace,v)),ST(map(rangespace,v)))
DiagonalInterlaceOperator(v::Vector,k...)=DiagonalInterlaceOperator(tuple(v...),k...)

for TYP in (:BandedOperator,:Operator)
    @eval Base.convert{T}(::Type{$TYP{T}},op::DiagonalInterlaceOperator)=
        DiagonalInterlaceOperator(map($TYP{T},op.ops),op.domainspace,op.rangespace)
end


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
