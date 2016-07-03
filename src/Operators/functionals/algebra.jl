## Operations
function *(A::Operator,b::Vector)
    @assert isafunctional(A)
    dotu(A[1:length(b)],b)
end
*(A::Operator,b::Fun) = promotedomainspace(A,space(b))*b.coefficients

*{T,D<:DefiniteIntegral,M<:Multiplication}(A::TimesFunctional{T,D,M},b::Fun) = bilinearform(A.op.f,b)
*{T,D<:DefiniteLineIntegral,M<:Multiplication}(A::TimesFunctional{T,D,M},b::Fun) = linebilinearform(A.op.f,b)


function *(c::Number,B::Operator)
    @assert isafunctional(B)
    c==1?B:ConstantTimesFunctional(c,B)
end
function *(B::Operator,c::Number)
    @assert isafunctional(B)
    c==1?B:ConstantTimesFunctional(c,B)
end
function /(B::Operator,c::Number)
    @assert isafunctional(B)
    c==1?B:ConstantTimesFunctional(1.0/c,B)
end
function *(B::Operator,O::TimesOperator)
    @assert isafunctional(B)
    TimesFunctional(B,O)  # Needed to avoid ambiguity
end
function *(B::Operator,O::BandedOperator)
    @assert isafunctional(B)
    if isconstop(O)
        promotedomainspace(B*convert(eltype(O),O),domainspace(O))
    else
        TimesFunctional(promotedomainspace(B,rangespace(O)),O)
    end
end

function -(B::Operator)
    @assert isafunctional(B)
    ConstantTimesFunctional(-1,B)
end


function -(A::Operator,B::Operator)
    @assert isafunctional(B)
    PlusFunctional([A,-B])
end
