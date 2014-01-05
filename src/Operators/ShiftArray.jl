export ShiftArray

type ShiftArray{T<:Number}
  data::Array{T,2}#TODO: probably should have more cols than rows
  colindex::Integer
end

ShiftArray(ind)=ShiftArray(Array(Float64,0,0),ind)

Base.setindex!(S::ShiftArray,x,k::Integer,j::Integer)=(S.data[k, j + S.colindex] = x)
Base.getindex(S::ShiftArray,k,j) = S.data[k, j + S.colindex]


function Base.resize!(S::ShiftArray,n::Integer,m::Integer)
    olddata = S.data
    S.data = zeros(n,m)
    
    if size(olddata,1) != 0
        S.data[1:size(olddata,1),1:m] = olddata[:,1:m]
    end
    
    S
end