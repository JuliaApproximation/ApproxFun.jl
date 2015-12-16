

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

function stride_addentries!(op,ri,ci,rs,cs,A,kr::UnitRange)
    r1=divrowrange(rs,ri,kr)
    if length(r1)>0
        addentries!(op,IndexStride(A,ri,ci,rs,cs),r1,:)
    end

    A
end



stride_addentries!(S::StrideOperator,A,kr::Range)=stride_addentries!(S.op,S.rowindex,S.colindex,S.rowstride,S.colstride,A,kr)


addentries!(S::StrideOperator,A,kr,::Colon)=stride_addentries!(S,A,kr)
domain(S::StrideOperator)=Any ##TODO: tensor product


## SliceOperator

# Some of this is verbatim from IndexSlice
immutable SliceOperator{T,B} <: BandedOperator{T}
    op::B
    rowindex::Int
    colindex::Int
    rowstride::Int
    colstride::Int

    function SliceOperator(o,r,c,rs,cs)
        @assert rs == cs
        @assert rs != 0
        @assert mod(r-c,rs)==0
        @assert mod(stride(o),rs)==0

        new(o,r,c,rs,cs)
    end
end

SliceOperator{T<:Number}(B::Operator{T},r,c,rs,cs)=SliceOperator{T,typeof(B)}(B,r,c,rs,cs)
SliceOperator{T<:Number}(B::Operator{T},r,c,rs)=SliceOperator{T,typeof(B)}(B,r,c,rs,rs)
SliceOperator{T<:Number}(B::Operator{T},r,c)=SliceOperator{T,typeof(B)}(B,r,c,1,1)


Base.convert{BT<:Operator}(::Type{BT},S::SliceOperator)=SliceOperator(convert(BandedOperator{eltype(BT)},S.op),
                                                                        S.rowindex,S.colindex,S.rowstride,S.colstride)

bandinds(S::SliceOperator)=(div(bandinds(S.op,1)+S.rowindex-S.colindex,S.rowstride),div(bandinds(S.op,2)+S.rowindex-S.colindex,S.rowstride))

function destride_addentries!(op,ri,ci,rs,cs,A,kr::UnitRange)
    r1=rs*kr[1]+ri:rs:rs*kr[end]+ri
    addentries!(op,IndexSlice(A,ri,ci,rs,cs),r1,:)
    A
end

function destride_addentries!(op,ri,ci,A,kr::UnitRange)
    r1=kr[1]+ri:kr[end]+ri
    addentries!(op,IndexSlice(A,ri,ci,1,1),r1,:)
    A
end

function destride_addentries!(S::SliceOperator,A,kr::Range)
    if S.rowstride==S.colstride==1
        destride_addentries!(S.op,S.rowindex,S.colindex,A,kr)
    else
        destride_addentries!(S.op,S.rowindex,S.colindex,S.rowstride,S.colstride,A,kr)
    end
end

addentries!(S::SliceOperator,A,kr,::Colon)=destride_addentries!(S,A,kr)
domain(S::SliceOperator)=domain(S.op)
domainspace(S::SliceOperator)=S.colindex==0&&S.colstride==1?domainspace(S.op):SliceSpace(domainspace(S.op),S.colindex,S.colstride)
rangespace(S::SliceOperator)=SliceSpace(rangespace(S.op),S.rowindex,S.rowstride)




## StrideFunctional


type StrideFunctional{T<:Number,B<:Functional} <: Functional{T}
    op::B
    rowindex::Int
    stride::Int
end

StrideFunctional{T<:Number}(B::Functional{T},r,rs)=StrideFunctional{T,typeof(B)}(B,r,rs)


Base.getindex{T<:Number}(op::StrideFunctional{T},kr::Range)=T[((k-op.rowindex)≥1 &&(k-op.rowindex)%op.stride==0)?op.op[fld(k-op.rowindex,op.stride)]:zero(T) for k=kr]

for TYP in (:Operator,:Functional)
    @eval Base.convert{T}(::Type{$TYP{T}},S::StrideFunctional)=StrideFunctional(convert(Functional{T},S.op),S.rowindex,S.stride)
end


##interlace block operators

iszerooperator(::ZeroOperator)=true
iszerooperator(::ZeroFunctional)=true
iszerooperator(A::ConstantOperator)=A.c==0.
iszerooperator(A)=false
function isboundaryrow(A,k)
    for j=1:size(A,2)
        if isa(A[k,j],Functional)
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
bandinds(M::InterlaceOperator,k::Integer)=k==1?(size(M.ops,k)*mapreduce(m->bandinds(m,k)-1,min,M.ops)+1):(size(M.ops,k)*mapreduce(m->bandinds(m,k)+1,max,M.ops)-1)
bandinds(M::InterlaceOperator)=bandinds(M,1),bandinds(M,2)

function addentries!(M::InterlaceOperator,A,kr::Range,::Colon)
    n=size(M.ops,1)
    for k=1:n,j=1:n
        stride_addentries!(M.ops[k,j],k-n,j-n,n,n,A,kr)
    end
    A
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
    if n==2 &&all(isconstop,A[1:end-1,1]) &&iscolop(A[end,1]) &&
            all(a->isa(a,Functional),A[1:end-1,2]) && isa(A[end,2],BandedOperator)
        return [Functional{TT}[BlockFunctional(convert(Number,A[k,1]),A[k,2]) for k=1:m-1]...;
                    BlockOperator(A[end,:])]
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


function addentries!(D::DiagonalInterlaceOperator,A,kr::Range,::Colon)
    n=length(D.ops)
    for k=1:n
        stride_addentries!(D.ops[k],k-n,k-n,n,n,A,kr)
    end
    A
end

domainspace(D::DiagonalInterlaceOperator)=D.domainspace
rangespace(D::DiagonalInterlaceOperator)=D.rangespace
