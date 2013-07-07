

export ShiftVector



## Interval

type ShiftVector{T<:Real}
    vector::Vector{T}
    index::Integer
end

ShiftVector(neg::Vector,nonneg::Vector)=ShiftVector([neg,nonneg],length(neg)+1)

Base.length(sl::ShiftVector)=length(sl.vector)
Base.getindex(sl::ShiftVector,k::Integer)=sl.vector[k+sl.index]
firstindex(sl::ShiftVector)=1-sl.index
lastindex(sl::ShiftVector)=length(sl)-sl.index

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