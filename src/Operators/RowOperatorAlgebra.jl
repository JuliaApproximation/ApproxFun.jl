
##PlusRowOperator

type PlusRowOperator{T<:Number,B<:RowOperator} <: RowOperator{T} 
    ops::Vector{B}
end

PlusRowOperator{B<:RowOperator}(ops::Vector{B})=PlusRowOperator{Float64,B}(ops)

Base.getindex(op::PlusRowOperator,k::Range1)=mapreduce(o->o[k],+,op.ops)

+(A::PlusRowOperator,B::PlusRowOperator)=PlusRowOperator([A.ops,B.ops])
+(A::PlusRowOperator,B::RowOperator)=PlusRowOperator([A.ops,B])
+(A::RowOperator,B::PlusRowOperator)=PlusRowOperator([A,B.ops])
+(A::RowOperator,B::RowOperator)=PlusRowOperator([A,B])



type ConstantTimesRowOperator{T<:Number,B<:RowOperator} <: RowOperator{T}
    c::T
    op::B
end

Base.getindex(op::ConstantTimesRowOperator,k::Range1)=op.c*op.op[k]

*(c::Number,B::RowOperator)=ConstantTimesRowOperator(c,B)
*(B::RowOperator,c::Number)=ConstantTimesRowOperator(c,B)
-{T<:Number}(B::RowOperator{T})=ConstantTimesRowOperator(-one(T),B)


-(A::RowOperator,B::RowOperator)=PlusRowOperator([A,-B])