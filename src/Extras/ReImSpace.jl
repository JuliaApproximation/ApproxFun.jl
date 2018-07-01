


## just takes realpart of operator
struct ReOperator{O,T} <: Operator{T}
    op::O
end

ReOperator(op)=ReOperator{typeof(op),Float64}(op)
convert(::Type{Operator{T}},R::ReOperator) where {T} = ReOperator{typeof(R.op),T}(R.op)

@wrapperstructure ReOperator
@wrapperspaces ReOperator





getindex(RI::ReOperator{O,T},k::Integer,j::Integer) where {O,T} =
    convert(T,real(RI.op[k,j]))

choosedomainspace(R::ReOperator,sp::Space) = choosedomainspace(R.op,sp)
for OP in (:promotedomainspace,:promoterangespace)
    @eval begin
        $OP(R::ReOperator,sp::UnsetSpace) = ReOperator($OP(R.op,sp))
        $OP(R::ReOperator,sp::Space) = ReOperator($OP(R.op,sp))
    end
end



# TODO: can't do this because UnsetSpace might change type
#real{T<:Real}(op::Operator{T})=op
real(op::Operator) = ReOperator(op)
