

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


divrowrange(rs,ri,r)=firstrw(rs,ri,r[1]):lastrw(rs,ri,r[end])

for op in (:firstrw,:lastrw,:divrowrange)
    @eval $op(S,k...)=$op(S.rowstride,S.rowindex,k...)
end


#S[rowstride*k + rowindex,colstride*j + colindex] == op[k,j]
#S[k,j] == A[k,j-k]
#A[rowstride*k + rowindex,colstride*j + colindex - k] == op[k,j]

function stride_addentries!(op,ri,ci,rs,cs,A,kr::UnitRange)
    r1=divrowrange(rs,ri,kr)
    
    addentries!(op,IndexStride(A,ri,ci,rs,cs),r1)
    
    A    
end



stride_addentries!(S::StrideOperator,A,kr::Range)=stride_addentries!(S.op,S.rowindex,S.colindex,S.rowstride,S.colstride,A,kr)


addentries!(S::StrideOperator,A,kr)=stride_addentries!(S,A,kr)
domain(S::StrideOperator)=Any ##TODO: tensor product


## DestrideOperator

# Some of this is verbatim from IndexDestride
immutable DestrideOperator{T<:Number,B<:Operator} <: BandedOperator{T}
    op::B
    rowindex::Int       
    colindex::Int       
    rowstride::Int
    colstride::Int
    
    function DestrideOperator(o,r,c,rs,cs)
        @assert rs == cs
        @assert rs != 0
        @assert mod(r-c,rs)==0        
        @assert mod(stride(o),rs)==0
        
        new(o,r,c,rs,cs)
    end
end

DestrideOperator{T<:Number}(B::Operator{T},r,c,rs,cs)=DestrideOperator{T,typeof(B)}(B,r,c,rs,cs)
DestrideOperator{T<:Number}(B::Operator{T},r,c,rs)=DestrideOperator{T,typeof(B)}(B,r,c,rs,rs)


bandinds(S::DestrideOperator)=(div(bandinds(S.op,1)+S.colindex-S.rowindex,S.rowstride),div(bandinds(S.op,2)+S.colindex-S.rowindex,S.rowstride))

function destride_addentries!(op,ri,ci,rs,cs,A,kr::UnitRange)
    r1=rs*kr[1]+ri:rs:rs*kr[end]+ri
    
    addentries!(op,IndexDestride(A,ri,ci,rs,cs),r1)
    
    A    
end

destride_addentries!(S::DestrideOperator,A,kr::Range)=destride_addentries!(S.op,S.rowindex,S.colindex,S.rowstride,S.colstride,A,kr)

addentries!(S::DestrideOperator,A,kr)=destride_addentries!(S,A,kr)
domain(S::DestrideOperator)=domain(S.op)
domainspace(S::DestrideOperator)=StrideSpace(domainspace(S.op),S.colindex,S.colstride)
rangespace(S::DestrideOperator)=StrideSpace(rangespace(S.op),S.rowindex,S.rowstride)




## StrideFunctional


type StrideFunctional{T<:Number,B<:Functional} <: Functional{T}
    op::B
    rowindex::Int
    stride::Int  
end

StrideFunctional{T<:Number}(B::Functional{T},r,rs)=StrideFunctional{T,typeof(B)}(B,r,rs)


Base.getindex{T<:Number}(op::StrideFunctional{T},kr::Range1)=[((k-op.rowindex)%op.stride==0)?op.op[fld(k-op.rowindex,op.stride)]:zero(T) for k=kr]



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
    else
        PiecewiseSpace(spl)
    end
end

function rangespace{T<:Operator}(A::Vector{T})
    @assert spacescompatible(A)
    
    spl=map(rangespace,A)
    if spacescompatible(spl)
        ArraySpace(first(spl),length(spl))
    else
        PiecewiseSpace(spl)
    end
end

function promotespaces{T<:Operator}(A::Array{T,2})
    A=copy(A)#TODO: promote might have different Array type
    for j=1:size(A,2)
        A[:,j]=promotedomainspace(A[:,j])
    end
    for k=1:size(A,1)
        A[k,:]=promoterangespace(vec(A[k,:]))
    end
    A
end


## Interlace operator

immutable InterlaceOperator{T} <: BandedOperator{T}
    ops::Matrix{BandedOperator{T}} 
    function InterlaceOperator(os)
        @assert size(os,1)==size(os,2)
        new(os)
    end
end
InterlaceOperator{T}(ops::Matrix{BandedOperator{T}})=InterlaceOperator{T}(ops)

bandinds(M::InterlaceOperator,k::Integer)=k==1?(size(M.ops,k)*mapreduce(m->bandinds(m,k),min,M.ops)):(size(M.ops,k)*mapreduce(m->bandinds(m,k),max,M.ops))
bandinds(M::InterlaceOperator)=bandinds(M,1),bandinds(M,2)

function addentries!(M::InterlaceOperator,A,kr::Range)
    n=size(M.ops,1)
    for k=1:n,j=1:n
        stride_addentries!(M.ops[k,j],k-n,j-n,n,n,A,kr) 
    end
    A
end



function interlace{T<:Operator}(A::Matrix{T})
    m,n=size(A)
    
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
    
    S=Array(Operator,br<m?br+1:br)
    
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

    for k=br+1:m
        Ap=vec(A[k,:])
        
        for j=1:n
            if !iszerooperator(A[k,j])  #not sure what promote does for constant operator
                op = StrideOperator(Ap[j],k-br-n,j-n,n)
            
                if !isdefined(S,br+1)
                    S[br+1] = op
                else
                    S[br+1] = S[br+1] + op
                end
            end            
        end
    end
    
    if br < m
        rsp=rangespace(A[br+1:end,1])    
        S[br+1]=SpaceOperator(S[br+1],dsp,rsp)
    end
    
    if(size(S,1) ==1)
        S[1]
    else
        S
    end
end


abstract AbstractDiagonalInterlaceOperator{T,B}<:BandedOperator{T}

immutable DiagonalInterlaceOperator{T<:Number,B<:Operator} <: AbstractDiagonalInterlaceOperator{T,B}
    ops::Vector{B}
end

DiagonalInterlaceOperator{B<:Operator}(v::Vector{B})=DiagonalInterlaceOperator{mapreduce(eltype,promote_type,v),B}(v)
DiagonalInterlaceOperator(v::Vector{Any})=DiagonalInterlaceOperator(Operator{mapreduce(eltype,promote_type,v)}[v...])



function bandinds(S::AbstractDiagonalInterlaceOperator)
    binds=map(bandinds,S.ops)
    bra=mapreduce(first,min,binds)
    brb=mapreduce(last,max,binds)    
    n=length(S.ops)
    n*bra,n*brb
end


function addentries!(D::AbstractDiagonalInterlaceOperator,A,kr::Range)
    n=length(D.ops)
    for k=1:n
        stride_addentries!(D.ops[k],k-n,k-n,n,n,A,kr)
    end
    A
end


for op in (:domainspace,:rangespace)
    @eval $op(D::DiagonalInterlaceOperator)=devec(map($op,D.ops))
end


