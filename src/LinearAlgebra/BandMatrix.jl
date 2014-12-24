

## 
# Represent a band matrix
# [ a_11 a_12 
#   a_21 a_22 a_23
#   a_31 a_32 a_33 a_34
#        a_42 a_43 a_44  ]
# ordering the data like
#       [ *     *       a_31    a_42
#         *      a_21   a_32    a_43 
#         a_11   a_22   a_33    A_44
#         a_12   a_23   a_44    *       ]
###


immutable BandMatrix{T}
    data::Matrix{T}
    l::Int # lower bandwidth ≥0
    u::Int # upper bandwidth ≥0
    m::Int #Number of columns
    function BandMatrix(data::Matrix{T},a,b,m)
        @assert size(data,1)==b+a+1
        new(data,a,b,m)
    end    
end


BandMatrix{T}(data::Matrix{T},a,b,m)=BandMatrix{T}(data,a,b,m)
BandMatrix{T}(::Type{T},a::Integer,b,n,m)=BandMatrix{T}(Array(T,b+a+1,n),a,b,m)
BandMatrix{T}(::Type{T},a::Integer,b,n)=BandMatrix{T}(Array(T,b+a+1,n),a,b,n)


for (op,bop) in ((:(Base.rand),:barand),(:(Base.zeros),:bazeros),(:(Base.ones),:baones))
    @eval begin
        $bop{T}(::Type{T},a::Integer,b,n,m)=BandMatrix($op(T,b+a+1,n),a,b,m)
        $bop{T}(::Type{T},a::Integer,b,n)=$bop(T,a,b,n,n)
        $bop(a::Integer,b,n,m)=$bop(Float64,a,b,n,m)
        $bop(a::Integer,b,n)=$bop(Float64,a,b,n,n)
    end
end


Base.size(A::BandMatrix,k)=ifelse(k==1,size(A.data,2),A.m)
Base.size(A::BandMatrix)=size(A.data,2),A.m
Base.eltype{T}(::BandMatrix{T})=T

usgetindex(A::BandMatrix,k,j::Integer)=A.data[j-k+A.l+1,k]
getindex{T}(A::BandMatrix{T},k::Integer,j::Integer)=(-A.l≤j-k≤A.u)?usgetindex(A,k,j):(j≤A.m?zero(T):throw(BoundsError()))
getindex(A::BandMatrix,kr::Range,j::Integer)=-A.l≤j-kr[end]≤j-kr[1]≤A.u?usgetindex(A,kr,j):[A[k,j] for k=kr]
getindex(A::BandMatrix,k::Integer,jr::Range)=[A[k,j] for j=jr]
getindex(A::BandMatrix,kr::Range,jr::Range)=[A[k,j] for k=kr,j=jr]
Base.full(A::BandMatrix)=A[1:size(A,1),1:size(A,2)]


# We turn off bound checking to allow nicer syntax without branching
#setindex!(A::BandMatrix,v,k::Integer,j::Integer)=((A.l≤j-k≤A.u)&&k≤A.n)?ussetindex!(A,v,k,j):throw(BoundsError())
#setindex!(A::BandMatrix,v,kr::Range,j::Integer)=(A.l≤j-kr[end]≤j-kr[1]≤A.u&&kr[end]≤A.n)?ussetindex!(A,v,kr,j):throw(BoundsError())


ibsetindex!(A::BandMatrix,v,k,j::Integer)=(@inbounds A.data[j-k+A.l+1,k]=v)
ibpluseq!(A::BandMatrix,v,k,j::Integer)=(@inbounds A.data[j-k+A.l+1,k]+=v)
setindex!(A::BandMatrix,v,k,j::Integer)=(A.data[j-k+A.l+1,k]=v)

function setindex!(A::BandMatrix,v,kr::Range,jr::Range)
    for j in jr
        A[kr,j]=slice(v,:,j)
    end
end
function setindex!(A::BandMatrix,v,k::Integer,jr::Range)
    for j in jr
        A[k,j]=v[j]
    end
end



function bamultiply!(C::BandMatrix,A::BandMatrix,B::BandMatrix)   
    n,m=size(C)
    for k=1:n  # ROWS
        for l=max(1,k-A.l):min(k+A.u,size(A,2)) # columns of A
            @inbounds Aj=A.data[l-k+A.l+1,k]
            
            shA=-l+B.l+1
            shB=-k+C.l+l-B.l
            @simd for j=max(1,k-C.l,l-B.l)+shA:min(B.u+l,n)+shA # columns of C/B
                @inbounds C.data[j+shB,k]+=Aj*B.data[j,l]
            end
        end
    end 
    C
end



function *{T,V}(A::BandMatrix{T},B::BandMatrix{V})
    if size(A,2)!=size(B,1)
        throw(DimensionMismatch("*"))
    end
    n=size(A,1)
    m=size(B,2)    
    bamultiply!(bazeros(promote_type(T,V),A.l+B.l,A.u+B.u,n,m),A,B)
end



