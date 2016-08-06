


## just takes realpart of operator
immutable ReOperator{O,T} <: Operator{T}
    op::O
end

ReOperator(op)=ReOperator{typeof(op),Float64}(op)
Base.convert{BT<:Operator}(::Type{BT},R::ReOperator)=ReOperator{typeof(R.op),eltype(BT)}(R.op)

for OP in (:rangespace,:domainspace,:bandinds)
    @eval $OP(R::ReOperator)=$OP(R.op)
end



getindex(RI::ReOperator,k::Integer,j::Integer) =
    real(RI.op[k,j])

choosedomainspace(R::ReOperator,sp::Space) = choosedomainspace(R.op,sp)
for OP in (:promotedomainspace,:promoterangespace)
    @eval begin
        $OP(R::ReOperator,sp::UnsetSpace) = ReOperator($OP(R.op,sp))
        $OP(R::ReOperator,sp::Space) = ReOperator($OP(R.op,sp))
    end
end



# TODO: can't do this because UnsetSpace might change type
#Base.real{T<:Real}(op::Operator{T})=op
Base.real(op::Operator)=ReOperator(op)
