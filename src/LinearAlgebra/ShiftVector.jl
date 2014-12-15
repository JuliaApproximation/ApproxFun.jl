

export ShiftVector,firstindex,lastindex,index



## TODO: Allow general T

type ShiftVector{T<:Number}
    vector::Vector{T}
    index::Int
end

ShiftVector(neg::Vector,nonneg::Vector)=ShiftVector([neg,nonneg],length(neg)+1)


for op = (:(Base.last),:(Base.first),:(Base.norm),:(Base.length))
    @eval ($op)(sv::ShiftVector) = ($op)(sv.vector)
end

Base.endof(sl::ShiftVector)=endof(sl.vector)-sl.index

index(sl::ShiftVector)=sl.index
firstindex(sl::ShiftVector)=1-sl.index
lastindex(sl::ShiftVector)=length(sl)-sl.index

range(sv::ShiftVector)=firstindex(sv):lastindex(sv)

Base.getindex(sl::ShiftVector,k::Integer)=sl.vector[k+sl.index]
Base.getindex(sl::ShiftVector,r::Range1)=sl.vector[r+sl.index]


Base.flipud(sl::ShiftVector)=ShiftVector(flipud(sl.vector),length(sl.vector)-sl.index+1)

##Assignment

Base.setindex!(sv::ShiftVector,x::Number,k::Integer)=sv.vector[k+sv.index]=x



for op = (:*,:.*,:./,:/)
    @eval ($op)(f::ShiftVector,c::Number) = ShiftVector(($op)(f.vector,c),f.index)
end 

-(f::ShiftVector)=ShiftVector(-f.vector,f.index)
-(c::Number,f::ShiftVector)=-(f-c)


for op = (:*,:.*,:+)
    @eval ($op)(c::Number,f::ShiftVector)=($op)(f,c)
end



##TODO: is being able to add vectors of different sizes desired behaviour?
for op = (:+,:-,:.*,:./)
    @eval begin
        function ($op)(f::ShiftVector,g::ShiftVector)
            @assert range(f) == range(g)
            @assert f.index == g.index
            
            ShiftVector(($op)(f.vector,g.vector),f.index)
        end
    end
end 



function Base.print(sl::ShiftVector)
    print("[")
    for k = firstindex(sl):-1
        print(sl[k],",")
    end
    print("_",sl[0],"_")
    for k = 1:lastindex(sl)
        print(",",sl[k])
    end    
    print("]")
end


#TODO f[0:end]

pad(f::ShiftVector,r::Range1)=ShiftVector(
    padleft(f[firstindex(f):-1],-first(r)),
    pad(f[0:lastindex(f)],last(r)+1))






##FFT That returns ShiftVector
function svfft(v::Vector)
        n=length(v)
        v=FFTW.fft(v)./n
        if mod(n,2) == 0
            ind=convert(Integer,n/2)
            v=alternatesign!(v)
            ShiftVector(v[ind+1:end],            
                        v[1:ind])            
        elseif mod(n,4)==3
            ind=convert(Integer,(n+1)/2)
            ShiftVector(-alternatesign!(v[ind+1:end]),     
                        alternatesign!(v[1:ind]))             
        else #mod(length(v),4)==1
            ind=convert(Integer,(n+1)/2)
            ShiftVector(alternatesign!(v[ind+1:end]),        
                        alternatesign!(v[1:ind]))             
        end
end





function isvfft(sv::ShiftVector)
        n=length(sv)
        ind = sv.index        
        
        #TODO: non normalized index ranges
        @assert n/2 <= ind <= n/2+1
        
        if mod(n,2) == 0
            v=alternatesign!([sv[0:lastindex(sv)],sv[firstindex(sv):-1]])
        elseif mod(n,4)==3
            v=[alternatesign!(sv[0:lastindex(sv)]),
                -alternatesign!(sv[firstindex(sv):-1])]    
        else #mod(length(v),4)==1
            v=[alternatesign!(sv[0:lastindex(sv)]),alternatesign!(sv[firstindex(sv):-1])]         
        end
        
        FFTW.ifft(n*v)
end



##interlace pos and neg coefficients
function interlace{T}(v::ShiftVector{T})
    fi=firstindex(v);li=lastindex(v)

    ret=zeros(T,(-fi>li)?-2*fi:2*li+1)
        
    
    j=1
    for k=0:li
        ret[j]=v[k]
        j+=2
    end
    
    j=2
    for k=-1:-1:fi
        ret[j]=v[k]
        j+=2
    end
    
    ret
end

function interlace{T<:Number}(v::Vector{Vector{T}})
    n=length(v)
    l=mapreduce(length,max,v)
    ret=zeros(T,n*l)
    
    for k=1:n
        ret[k:n:k+(length(v[k])-1)*n]=v[k]
    end
    ret
end

function interlace(v::Vector{Any})
    #determine type
    T=Float64
    for vk in v
        if isa(vk,Vector{Complex{Float64}})
            T=Complex{Float64}
        end
    end
    b=Array(Vector{T},length(v))
    for k=1:length(v)
        b[k]=v[k]
    end
    interlace(b)
end

interlace{T}(a::Vector{T},b::Vector{T})=interlace(Vector{T}[a,b])

deinterlace(v::Vector)=ShiftVector(flipud(v[2:2:end]),v[1:2:end])



function Base.resize!(c::ShiftVector,k::Range1)
   fi=firstindex(c)

   splice!(c.vector,last(k)+c.index+1:length(c.vector)) 
   splice!(c.vector,1:first(k)+c.index-1)
   c.index += fi - first(k) 
   
   c
end


## chop
##TODO: what if all entries are zero
function chop!(c::ShiftVector,tol::Real)
    @assert tol > 0

    for k=range(c)
        if abs(c[k]) > tol
            resize!(c,k:lastindex(c))
            break
        end
    end
    
    for k=lastindex(c):-1:firstindex(c)
        if abs(c[k]) > tol
            resize!(c,firstindex(c):k)
            break
        end
    end    
    
    c
end

Base.chop(c::ShiftVector,tol::Real)=chop!(deepcopy(c),tol)
