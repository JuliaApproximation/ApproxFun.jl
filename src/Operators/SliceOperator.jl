

export SliceOperator


type SliceOperator{T<:Number,B<:Operator{T}} <: BandedOperator{T}
    op::B
    rowindex::Int       
    colindex::Int       
    rowstride::Int
    colstride::Int
    
    function SliceOperator(o,r,c,rs,cs)
        @assert rs == 1 ##TODO: general steps
        @assert cs == 1        
        
        new(o,r,c,rs,cs)
    end
end



SliceOperator{T<:Number}(B::Operator{T},r,c,rs,cs)=SliceOperator{T,typeof(B)}(B,r,c,rs,cs)
SliceOperator{T<:Number}(B::Operator{T},r,c,rs)=SliceOperator{T,typeof(B)}(B,r,c,rs,rs)

function bandinds(S::SliceOperator)
    br=bandinds(S.op)
    
    sh=S.colindex-S.rowindex
    
    (min(br[1]-sh,0),max(0,br[end]-sh))
end
 
function addentries!(S::SliceOperator,A,kr)
    sh=S.rowindex-S.colindex
    Ash=ShiftArray(A.data,A.rowindex-S.rowindex+1,A.colindex+sh)
    addentries!(S.op,Ash,kr+S.rowindex-1)
    A
end


##TODO Slice spaces
domain(S::SliceOperator)=AnySpace()


Base.ndims(::BandedOperator)=2
function Base.getindex(B::BandedOperator,k::Base.FloatRange,j::Base.FloatRange)
    @assert last(k)==Inf && last(j)==Inf
    
    SliceOperator(B,convert(Integer,first(k)),convert(Integer,first(j)),1,1)
end