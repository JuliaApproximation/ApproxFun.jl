

export StrideOperator,StrideFunctional



#S[rowstride*k + rowindex,colstride*j + colindex] == op[k,j]
#S[k,j] == op[(k-rowindex)/rowstride,(j-colindex)/colstride]
type StrideOperator{T<:Number,B<:Operator} <: BandedOperator{T}
    op::B
    rowindex::Int       
    colindex::Int       
    rowstride::Int
    colstride::Int
    
    function StrideOperator(o,r,c,rs,cs)
        @assert abs(rs) == abs(cs)
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

function stride_pospos_addentries!(S::StrideOperator,A::ShiftArray,kr::Range1)
    r1=divrowrange(S,kr)

    B1=BandedArray(S.op,r1)
    B=BandedArray(A)
    
    for k=r1, j=columnrange(B1.data)+k
        B[S.rowstride*k + S.rowindex,S.colstride*j + S.colindex] += B1.data[k,j-k]
    end
    
    A
end

addentries!(S::StrideOperator,A,kr)=stride_pospos_addentries!(S,A,kr)
domain(S::StrideOperator)=Any ##TODO: tensor product


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
        VectorDomainSpace(first(spl),length(spl))
    else
        PiecewiseSpace(spl)
    end
end

function rangespace{T<:Operator}(A::Array{T})
    @assert spacescompatible(A)
    
    spl=map(rangespace,A[:,1])
    if spacescompatible(spl)
        VectorDomainSpace(first(spl),length(spl))
    else
        PiecewiseSpace(spl)
    end
end

function promotespaces{T<:Operator}(A::Array{T,2})
    A=copy(A)
    for j=1:size(A,2)
        A[:,j]=promotedomainspace(A[:,j])
    end
    for k=1:size(A,1)
        A[k,:]=promoterangespace(vec(A[k,:]))
    end
    A
end

function interlace{T<:Operator}(A::Array{T,2})
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


