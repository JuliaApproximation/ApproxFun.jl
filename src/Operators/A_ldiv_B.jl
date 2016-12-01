## Linear Solve
doc"""
    \(A,b;tolerance=tol,maxlength=n)

solves a linear equation, usually differential equation, where `A` is an operator
or array of operators and `b` is a `Fun` or array of funs.  The result `u`
will approximately satisfy `A*u = b`.
"""
function \(A::Operator,b;kwds...)
    if isambiguous(domainspace(A))
        A=choosespaces(A,b)
        if isambiguous(domainspace(A))
            error("Cannot infer spaces")
        end
        \(A,b;kwds...)
    else
        \(qrfact(A),b;kwds...)
    end
end

A_ldiv_B_coefficients(A::Operator,b;kwds...) = A_ldiv_B_coefficients(qrfact(A),b;kwds...)


#TODO: Remove these when interlace is automatic
for TYP in (:Vector,:Matrix)
    @eval begin
        \{OO<:Operator,T}(A::$TYP{OO},b::AbstractVecOrMat{T};kwds...) =
            \(interlace(A),b;kwds...)
        \{OO<:Operator}(A::$TYP{OO},b;kwds...) = \(interlace(A),b;kwds...)
    end
end
A_ldiv_B_coefficients{OO<:Operator}(A::Array{OO},b;kwds...) =
    A_ldiv_B_coefficients(interlace(A),b;kwds...)
