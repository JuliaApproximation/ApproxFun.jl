export OperatorFunction


@compat abstract type OperatorFunction{BT,FF,T} <: Operator{T} end

immutable ConcreteOperatorFunction{BT<:Operator,FF,T} <: OperatorFunction{BT,FF,T}
    op::BT
    f::FF
end

ConcreteOperatorFunction(op::Operator,f::Function) =
    ConcreteOperatorFunction{typeof(op),typeof(f),eltype(op)}(op,f)
OperatorFunction(op::Operator,f::Function) = ConcreteOperatorFunction(op,f)

for op in (:domainspace,:rangespace,:domain,:bandinds)
    @eval begin
        $op(OF::ConcreteOperatorFunction) = $op(OF.op)
    end
end

function getindex(OF::ConcreteOperatorFunction,k::Integer,j::Integer)
    @assert isdiag(OF.op)
    if k==j
        OF.f(OF.op[k,k])::eltype(OF)
    else
        zero(eltype(OF))
    end
end

function Base.convert{T}(::Type{Operator{T}},D::ConcreteOperatorFunction)
    if T==eltype(D)
        D
    else
        ConcreteOperatorFunction{typeof(D.op),T}(D.op,D.f)
    end
end


for OP in (:(Base.inv),:(Base.sqrt))
    @eval begin
        $OP(D::DiagonalOperator) = OperatorFunction(D,$OP)
        $OP(C::ConstantTimesOperator) = $OP(C.λ)*$OP(C.op)
        function $OP(D::ConcreteOperatorFunction)
            @assert isdiag(D)
            OperatorFunction(D.op,x->$OP(D.f(x)))
        end
        $OP(A::Operator) = isdiag(A) ? OperatorFunction(A,$OP) :
                                error("Not implemented.")
    end
end


Base.sqrt(S::SpaceOperator) = SpaceOperator(sqrt(S.op),S.domainspace,S.rangespace)
Base.inv(S::SpaceOperator) = SpaceOperator(inv(S.op),S.rangespace,S.domainspace)
