## Linear Solve
doc"""
    `linsolve(A,b;tolerance=tol,maxlength=n)`

solves a linear equation, usually differential equation, where `A` is an operator
or array of operators and `b` is a `Fun` or array of funs.  The result `u`
will approximately satisfy `A*u = b`.
"""
function linsolve(A::Operator,b;kwds...)
    if isambiguous(domainspace(A))
        A=choosespaces(A,b)
        if isambiguous(domainspace(A))
            error("Cannot infer spaces")
        end
        linsolve(A,b;kwds...)
    else
        linsolve(qrfact(A),b;kwds...)
    end
end

linsolve{OO<:Operator}(A::Array{OO},b;kwds...) = linsolve(interlace(A),b;kwds...)


for p in (1,2)
    @eval begin
        \{T<:Operator}(A::Array{T,$p},b::Union{Array,Number,Fun}) = linsolve(A,b)
        \{T<:Operator,F<:Fun}(A::Array{T,$p},b::Vector{F}) = linsolve(A,b)
        \{T<:Operator,F<:Number}(A::Array{T,$p},b::Vector{F}) = linsolve(A,b)
    end
end
\(A::Operator,b) = linsolve(A,b)
