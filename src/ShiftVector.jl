

export ShiftVector



## Interval

type ShiftVector{T<:Number}
    vector::Vector{T}
    index::Integer
end

ShiftVector(neg::Vector,nonneg::Vector)=ShiftVector([neg,nonneg],length(neg)+1)

Base.length(sl::ShiftVector)=length(sl.vector)
Base.getindex(sl::ShiftVector,k::Integer)=sl.vector[k+sl.index]
firstindex(sl::ShiftVector)=1-sl.index
lastindex(sl::ShiftVector)=length(sl)-sl.index


for op = (:*,:.*,:./,:/)
    @eval ($op)(f::ShiftVector,c::Number) = ShiftVector(($op)(f.vector,c),f.index)
end 

-(f::ShiftVector)=ShiftVector(-f.vector,f.index)
-(c::Number,f::ShiftVector)=-(f-c)


for op = (:*,:.*,:+)
    @eval ($op)(c::Number,f::ShiftVector)=($op)(f,c)
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



##FFT That returns ShiftVector
function svfft(v)
        n=length(v)
        v=FFTW.fft(v)./n
        if(mod(length(v),2) == 0)
            ind=convert(Integer,length(v)/2)
            v=(-1)^n*alternatingvector(n).*v          
            ShiftVector(v[ind+1:end],            
                        v[1:ind])            
        elseif(mod(length(v),4)==3)
            ind=convert(Integer,(length(v)+1)/2)
            ShiftVector(-alternatingvector(n-ind).*v[ind+1:end],            
                        alternatingvector(ind).*v[1:ind])                
        else #mod(length(v),4)==1
            ind=convert(Integer,(length(v)+1)/2)
            ShiftVector(alternatingvector(n-ind).*v[ind+1:end],            
                        alternatingvector(ind).*v[1:ind])             
        end
end
