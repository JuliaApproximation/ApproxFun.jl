

export MutableAlmostBandedOperator, adaptiveqr!



function indexrange(b::BandedBelowOperator,k::Integer)
    ret = bandrange(b) + k
  
    (ret[1] < 1) ? (1:ret[end]) : ret
end

index(b::BandedBelowOperator)=1-bandrange(b)[1]



## MutableAlmostBandedOperator


type MutableAlmostBandedOperator{T<:Number,M<:BandedOperator} <: BandedBelowOperator
    bc::RowOperator
    op::M
    data::ShiftArray{T}  
    filldata::Vector{T}
    
    bcdata::Vector
    bcfilldata::T
    
    datalength::Integer
end


#TODO: index(op) + 1 -> length(bc) + index(op)
function MutableAlmostBandedOperator(bc::RowOperator,op::BandedOperator)
    data = ShiftArray(index(op)+1)
    bcdata = bc[1:bandrange(op)[end]+1]        
                
    MutableAlmostBandedOperator(bc,op,data,Float64[],bcdata,1.,0 )
end

function MutableAlmostBandedOperator{T<:Operator}(B::Vector{T})
    @assert length(B) == 2 ##TODO: make general
    @assert typeof(B[1]) <: RowOperator
    @assert typeof(B[2]) <: BandedOperator
    
    MutableAlmostBandedOperator(B[1],B[2])
end


index(B::MutableAlmostBandedOperator)=index(B.op)
numbcs(B::MutableAlmostBandedOperator)=1

# for bandrange, we save room for changed entries during Givens
bandrange(B::MutableAlmostBandedOperator)=(bandrange(B.op)[1]-1):(bandrange(B.op)[end]-bandrange(B.op)[1])
datalength(b::MutableAlmostBandedOperator)=b.datalength




function Base.getindex(B::MutableAlmostBandedOperator,kr::Range1,jr::Range1)
    ret = spzeros(length(kr),length(jr))
    
    
    for k = kr
        if k <= datalength(B) || k == 1 #bc
            for j=jr
                ret[k,j] = B[k,j]
            end
        else
            ir = bandrange(B) + k
            
            for j=max(ir[1],jr[1]):min(ir[end],jr[end])
                ret[k-kr[1]+1,j-jr[1]+1] = B[k,j]
            end
        end
    end
    
    ret
end


function Base.getindex(B::MutableAlmostBandedOperator,k::Integer,j::Integer)  
    ir = bandrange(B) + k
    nbc = numbcs(B)
    
    if k <= nbc
        (j <= length(B.bcdata))? B.bcdata[j] : B.bcfilldata*B.bc[j]
    elseif k-nbc <= datalength(B) && j <= ir[end] && ir[1] <= j
        B.data[k-nbc,j-k]
    elseif k-nbc <= datalength(B) && j > ir[end]
        B.filldata[k-nbc]*B.bc[j]
    else
        B.op[k-nbc,j]##TODO: Slow
    end
end


getindex!(b::MutableAlmostBandedOperator,kr::Range1,jr::Range1)=resizedata!(b,kr[end])[kr,jr]
getindex!(b::MutableAlmostBandedOperator,kr::Integer,jr::Integer)=resizedata!(b,kr)[kr,jr]

function resizedata!{T<:Number,M<:BandedOperator}(B::MutableAlmostBandedOperator{T,M},n::Integer)
    l = datalength(B)
    if n > l
        resize!(B.data,2n,length(bandrange(B)))
        resize!(B.filldata,2n)  
        for j=l+1:n
            B.filldata[j]=0.
        end
        
        addentries!(B.op,B.data,l+1:n)
        
        B.datalength = n
    end
    
    B
end

function Base.setindex!(B::MutableAlmostBandedOperator,x,k::Integer,j::Integer)
    if k==1
        B.bcdata[k,j] = x
    else
        resizedata!(B,k)      
        B.data[k-numbcs(B),j-k] = x
    end
    x
end

function setfilldata!(B::MutableAlmostBandedOperator,x,k::Integer,j::Integer)
    @assert j == 1 ##TODO: multiple bcs
    if k==1
        B.bcfilldata = x
    else
        B.filldata[k-numbcs(B)] = x
    end
end

getfilldata(B::MutableAlmostBandedOperator,k::Integer,j::Integer)=(k==1)?B.bcfilldata:B.filldata[k-numbcs(B)]



##TODO: decide adaptive resize
function givensreduce!(B::MutableAlmostBandedOperator,v::Vector,k1::Integer,k2::Integer,j1::Integer)
    a=B[k1,j1];b=B[k2,j1]
    sq=sqrt(a*a + b*b)
    a=a/sq;b=b/sq
    
    v[k1],v[k2] = a*v[k1] + b*v[k2],-b*v[k1] + a*v[k2]    
    
    
    #TODO: Assuming that left rows are already zero
    
    for j = j1:indexrange(B,k1)[end]
        B1 = B[k1,j]
        B2 = B[k2,j]
        
        B[k1,j],B[k2,j]= a*B1 + b*B2,-b*B1 + a*B2
    end
    
    for j=indexrange(B,k1)[end]+1:indexrange(B,k2)[end]
        B1 = B[k1,j]
        B2 = B[k2,j]
        
        B[k2,j]=a*B2 - b*B1
    end
    
    B1 = getfilldata(B,k1,1)
    B2 = getfilldata(B,k2,1)
    
    setfilldata!(B, a*B1 + b*B2,k1,1)
    setfilldata!(B,-b*B1 + a*B2,k2,1)    
    
    
    #TODO: assert that the remaining of k2 are zero
    B
end

function givensreduce!(B::MutableAlmostBandedOperator,v::Vector,k1::Range1,j1::Integer)
  for k=k1[2]:k1[end]
    givensreduce!(B,v,k1[1],k,j1)
  end
end

givensreduce!(B::MutableAlmostBandedOperator,v::Vector,j::Integer)=givensreduce!(B,v,j:(j-bandrange(B)[1]),j)


function backsubstitution!(B::MutableAlmostBandedOperator,u)
    n=length(u)
    b=bandrange(B)[end]
    
    
    
    for k=n:-1:n-b
        for j=k+1:n
            u[k]-=B[k,j]*u[j]
        end
          
        u[k] /= B[k,k]
    end
    
    pk = 0
    for k=n-b-1:-1:1
        pk = u[k+b+1]*B.bc[k+b+1] + pk
        for j=k+1:k+b
            u[k]-=B[k,j]*u[j]
        end
        
        u[k] -= getfilldata(B,k,1)*pk
          
        u[k] /= B[k,k]
    end
  
  u
end



adaptiveqr{T<:Operator}(B::Vector{T},v::Vector) = adaptiveqr!(MutableAlmostBandedOperator(B),copy(v))
function adaptiveqr!(B::MutableAlmostBandedOperator,v::Vector)
    u=[v,zeros(100)]
    
    l = length(v) + 100  
    resizedata!(B,l)
    
    
    j=1
    b=-bandrange(B)[1]
    
    tol=.001eps()
    
    
    ##TODO: check js
    while norm(u[j:j+b-1]) > tol
        if j + b == l
            u = [u,zeros(l)]      
            l *= 2
            resizedata!(B,l)
        end
        
        
        givensreduce!(B,u,j)
        j+=1
    end
  
    backsubstitution!(B,u[1:j-1])
end


