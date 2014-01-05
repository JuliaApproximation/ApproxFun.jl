

function Base.exp(f::IFun)
    xm=indmax(f)
    B=EvaluationOperator(xm,f.domain)
    A=[B,DifferentialOperator([-diff(f),1.])]
    A\[exp(f[xm]),0.]    
end