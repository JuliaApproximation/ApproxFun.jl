## Linear Solve
doc"""
    \(A,b;tolerance=tol,maxlength=n)

solves a linear equation, usually differential equation, where `A` is an operator
or array of operators and `b` is a `Fun` or array of funs.  The result `u`
will approximately satisfy `A*u = b`.
"""
function \(A::Operator,b::Fun;kwds...)
    if isambiguous(domainspace(A))
        A=choosespaces(A,b)
        if isambiguous(domainspace(A))
            error("Cannot infer spaces")
        end
        \(A,b;kwds...)
    else
        Fun(domainspace(A),
            A_ldiv_B_coefficients(A,coefficients(b,rangespace(A));kwds...))
    end
end

\(A::Operator,b::StridedVector;kwds...) = \(A,Fun(b,rangespace(A));kwds...)
\(A::Operator,b::AbstractVector;kwds...) = \(A,Fun(b,rangespace(A));kwds...)
\(A::Operator,b;kwds...) = \(A,Fun(b,rangespace(A));kwds...)

# Solve each column separately
function \(A::Operator,B::AbstractMatrix;kwds...)
    ds=domainspace(A)
    ret=Array(Fun{typeof(ds),promote_type(mapreduce(eltype,promote_type,B),eltype(ds))},
              1,size(B,2))
    for j=1:size(B,2)
        ret[:,j]=\(A,B[:,j];kwds...)
    end
    demat(ret)
end

A_ldiv_B_coefficients(A::Operator,b;kwds...) = A_ldiv_B_coefficients(qrfact(A),b;kwds...)


#TODO: Remove these when interlace is automatic
for TYP in (:Vector,:Matrix)
    @eval begin
        \{OO<:Operator}(A::$TYP{OO},b::StridedVecOrMat;kwds...) =
            \(interlace(A),b;kwds...)
        \{OO<:Operator}(A::$TYP{OO},b::AbstractVecOrMat;kwds...) =
            \(interlace(A),b;kwds...)
        \{OO<:Operator}(A::$TYP{OO},b::Fun;kwds...) =
            \(interlace(A),b;kwds...)
    end
end
A_ldiv_B_coefficients{OO<:Operator}(A::Array{OO},b;kwds...) =
    A_ldiv_B_coefficients(interlace(A),b;kwds...)
