
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



type FunctionalTimesBandedOperator{T<:Number,A<:Functional,B<:BandedOperator} <: Functional{T}
    functional::A
    op::B
end



FunctionalTimesBandedOperator{T<:Number}(A::Functional{T},B::BandedOperator{T})=FunctionalTimesBandedOperator{T,typeof(A),typeof(B)}(A,B)


function Base.getindex{T<:Number}(f::FunctionalTimesBandedOperator{T},jr::Range)#j is columns
    bi=ApproxFun.bandinds(f.op)
    B=BandedArray(f.op,(jr[1]-bi[end]):(jr[end]-bi[1]))
    r=zeros(T,length(jr))
    for j in jr, k=j-bi[end]:j-bi[1]
        if k>=1
            r[j-jr[1]+1]+=f.functional[k]*B[k,j]
        end
    end
    r
end


## Operations
*(A::Functional,b::Vector)=dot(A[1:length(b)],b)
*(A::Functional,b::IFun)=A*b.coefficients


*(c::Number,B::Functional)=ConstantTimesFunctional(c,B)
*(B::Functional,c::Number)=ConstantTimesFunctional(c,B)
*(B::Functional,O::TimesOperator)=FunctionalTimesBandedOperator(B,O)
*(B::Functional,O::BandedOperator)=FunctionalTimesBandedOperator(B,O)

-{T<:Number}(B::Functional{T})=ConstantTimesFunctional(-one(T),B)


-(A::Functional,B::Functional)=PlusFunctional([A,-B])