## Linear Solve


linsolve(A::Operator,b;kwds...) = linsolve(qrfact(A),b;kwds...)
linsolve{OO<:Operator}(A::Array{OO},b;kwds...) = linsolve(qrfact(A),b;kwds...)

for p in (1,2)
    @eval begin
        \{T<:Operator}(A::Array{T,$p},b::Union{Array,Number,Fun}) = linsolve(A,b)
        \{T<:Operator,F<:Fun}(A::Array{T,$p},b::Vector{F}) = linsolve(A,b)
        \{T<:Operator,F<:Number}(A::Array{T,$p},b::Vector{F}) = linsolve(A,b)
    end
end
\(A::Operator,b) = linsolve(A,b)
