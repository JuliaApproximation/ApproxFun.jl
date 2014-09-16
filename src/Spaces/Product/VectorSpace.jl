immutable VectorSpace{d} <: FunctionSpace
end

typealias ScalarSpace VectorSpace{1}

=={d}(::VectorSpace{d},::VectorSpace{d})=true