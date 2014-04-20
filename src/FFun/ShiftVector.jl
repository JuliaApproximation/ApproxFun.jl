

export ShiftVector,firstindex,lastindex,index



## Interval

type ShiftVector{T<:Number}
    vector::Vector{T}
    index::Integer
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


Base.flipud(sl::ShiftVector)=ShiftVector(flipud(sl[1:end]),flipud(sl[firstindex(sl):0]))

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
            v=alternatingvector(n).*v          
            ShiftVector(v[ind+1:end],            
                        v[1:ind])            
        elseif mod(n,4)==3
            ind=convert(Integer,(n+1)/2)
            ShiftVector(-alternatingvector(n-ind).*v[ind+1:end],            
                        alternatingvector(ind).*v[1:ind])                
        else #mod(length(v),4)==1
            ind=convert(Integer,(n+1)/2)
            ShiftVector(alternatingvector(n-ind).*v[ind+1:end],            
                        alternatingvector(ind).*v[1:ind])             
        end
end





function isvfft(sv::ShiftVector)
        n=length(sv)
        ind = sv.index        
        
        #TODO: non normalized index ranges
        @assert n/2 <= ind <= n/2+1
        
        if mod(n,2) == 0
            v=alternatingvector(n).*[sv[0:lastindex(sv)],sv[firstindex(sv):-1]]      
        elseif mod(n,4)==3
            v=[alternatingvector(n-ind+1).*sv[0:lastindex(sv)],
                -alternatingvector(ind-1).*sv[firstindex(sv):-1]]    
        else #mod(length(v),4)==1
            v=[alternatingvector(n-ind+1).*sv[0:lastindex(sv)],alternatingvector(ind-1).*sv[firstindex(sv):-1]]         
        end
        
        FFTW.ifft(n*v)
end



##interlace pos and neg coefficients
function interlace{T}(v::ShiftVector{T})
    ret=Array(T,length(v.vector))
    ret[1:2:end]=v.vector[v.index:end]
    ret[2:2:end]=v.vector[v.index-1:-1:1]
    ret
end

deinterlace(v::Vector)=ShiftVector(flipud(v[2:2:end]),v[1:2:end])

