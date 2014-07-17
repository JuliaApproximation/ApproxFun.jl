
##PlusFunctional

type PlusFunctional{T<:Number,B<:Functional} <: Functional{T} 
    ops::Vector{B}
end

PlusFunctional{B<:Functional}(ops::Vector{B})=PlusFunctional{Float64,B}(ops)

Base.getindex(op::PlusFunctional,k::Range1)=mapreduce(o->o[k],+,op.ops)

+(A::PlusFunctional,B::PlusFunctional)=PlusFunctional([A.ops,B.ops])
+(A::PlusFunctional,B::Functional)=PlusFunctional([A.ops,B])
+(A::Functional,B::PlusFunctional)=PlusFunctional([A,B.ops])
+(A::Functional,B::Functional)=PlusFunctional([A,B])



type ConstantTimesFunctional{T<:Number,B<:Functional} <: Functional{T}
    c::T
    op::B
end

Base.getindex(op::ConstantTimesFunctional,k::Range1)=op.c*op.op[k]


## Operations
*(A::Functional,b::Vector)=dot(A[1:length(b)],b)
*(A::Functional,b::IFun)=A*b.coefficients


*(c::Number,B::Functional)=ConstantTimesFunctional(c,B)
*(B::Functional,c::Number)=ConstantTimesFunctional(c,B)
-{T<:Number}(B::Functional{T})=ConstantTimesFunctional(-one(T),B)


-(A::Functional,B::Functional)=PlusFunctional([A,-B])