## Operations
*(A::Functional,b::Vector)=dotu(A[1:length(b)],b)
*(A::Functional,b::Fun)=promotedomainspace(A,space(b))*b.coefficients

*{T,D<:DefiniteIntegral,M<:AbstractMultiplication}(A::TimesFunctional{T,D,M},b::Fun) = dotu(A.op.f,b)
*{T,D<:DefiniteLineIntegral,M<:AbstractMultiplication}(A::TimesFunctional{T,D,M},b::Fun) = linedotu(A.op.f,b)


*(c::Number,B::Functional)=c==1?B:ConstantTimesFunctional(c,B)
*(B::Functional,c::Number)=c==1?B:ConstantTimesFunctional(c,B)
/(B::Functional,c::Number)=c==1?B:ConstantTimesFunctional(1.0/c,B)
*(B::Functional,O::TimesOperator)=TimesFunctional(B,O)  # Needed to avoid ambiguity
function *(B::Functional,O::BandedOperator)
    if isconstop(O)
        promotedomainspace(B*convert(eltype(O),O),domainspace(O))
    else
        TimesFunctional(promotedomainspace(B,rangespace(O)),O)
    end
end

-(B::Functional)=ConstantTimesFunctional(-1,B)


-(A::Functional,B::Functional)=PlusFunctional([A,-B])
