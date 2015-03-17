

export ShiftVector,firstindex,lastindex,index



## TODO: Allow general T

type ShiftVector{T<:Number}
    vector::Vector{T}
    index::Int
end

ShiftVector(neg::Vector,nonneg::Vector)=ShiftVector([neg;nonneg],length(neg)+1)


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


Base.flipud(sl::ShiftVector)=ShiftVector(flipdim(sl.vector,1),length(sl.vector)-sl.index+1)

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



deinterlace(v::Vector)=ShiftVector(flipdim(v[2:2:end],1),v[1:2:end])



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
