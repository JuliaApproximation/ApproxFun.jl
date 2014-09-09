immutable VectorSpace <: OperatorSpace
    dimension::Int
end

==(a::VectorSpace,b::VectorSpace)= a.dimension==b.dimension