# Operators that are diagonal in a dimension can be solved in O(n) operations
# isdiagop is used to inspect an operator to see if it is diagonal


Base.isdiag{T<:Number}(B::Operator{T}) = bandinds(B)==(0,0)


# multivariate case
isdiagop(K::KroneckerOperator,k) = isdiag(K.ops[k])
isdiagop(S::WrapperOperator,k) = isdiagop(S.op,k)
isdiagop(A::Union{PlusOperator,TimesOperator},k) = all(op->isdiagop(op,k),A.ops)


# diagop gets out the op corresponding to the k-th column
# TODO: diag in x direction


diagop(A::KroneckerOperator,col)=A.ops[2]


diagop(A::PlusOperator,col) = mapreduce(op->diagop(op,col),+,A.ops)
diagop(A::TimesOperator,col) = mapreduce(op->diagop(op,col),*,A.ops)
diagop(A::SpaceOperator,col) = diagop(A.op,col)
# TODO: Space for ConstantOperator
diagop(A::ConstantOperator,col) = ConstantOperator(A.c)
diagop(A::ConstantTimesOperator,col) = A.c*diagop(A.op,col)




diagop(L,k)=error("Override diagop for "*string(typeof(L)))
