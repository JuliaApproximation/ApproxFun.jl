## Linear Solve
for TYP in (:Fun,:StridedVector,:AbstractVector,:Any)
    @eval function \(A::Operator,b::$TYP;kwds...)
        if isambiguous(domainspace(A))
            A=choosespaces(A,b)
            if isambiguous(domainspace(A))
                error("Cannot infer spaces")
            end
            \(A,b;kwds...)
        else
            Fun(domainspace(A),
                ldiv_coefficients(A,coefficients(b,rangespace(A));kwds...))
        end
    end
end

"""
    \\(A::Operator,b;tolerance=tol,maxlength=n)

solves a linear equation, usually differential equation, where `A` is an operator
or array of operators and `b` is a `Fun` or array of funs.  The result `u`
will approximately satisfy `A*u = b`.
"""
\(::Operator,_)

# Solve each column separately
function \(A::Operator,B::AbstractMatrix;kwds...)
    ds=domainspace(A)
    if isambiguous(ds)
        return choosespaces(A,B[:,1])\B
    end

    ret=Matrix{VFun{typeof(ds),
               promote_type(eltype(A),mapreduce(eltype,promote_type,B))}}(1,size(B,2))

    QR = factorize(A) # reuse computation
    for j=1:size(B,2)
        ret[:,j]=\(QR,B[:,j];kwds...)
    end
    Fun(ret)
end

\(A::Operator,B::MatrixFun;kwds...) = \(A,Array(B);kwds...)

ldiv_coefficients(A::Operator,b;kwds...) = ldiv_coefficients(qrfact(A),b;kwds...)

\(A::Operator,B::Operator) = TimesOperator(inv(A),B)

#TODO: Remove these when interlace is automatic
for TYP in (:Vector,:Matrix)
    @eval begin
        \(A::$TYP{OO},b::StridedVecOrMat;kwds...) where {OO<:Operator} =
            \(interlace(A),b;kwds...)
        \(A::$TYP{OO},b::AbstractVecOrMat;kwds...) where {OO<:Operator} =
            \(interlace(A),b;kwds...)
        \(A::$TYP{OO},b::Fun;kwds...) where {OO<:Operator} =
            \(interlace(A),b;kwds...)
    end
end
ldiv_coefficients(A::AbstractArray{OO},b;kwds...) where {OO<:Operator} =
    ldiv_coefficients(interlace(A),b;kwds...)
