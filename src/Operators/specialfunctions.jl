

function Base.exp(f::IFun)
    xm=indmax(f)
    B=EvaluationOperator(f.domain,xm)
    D=diff(f.domain)
    A=[B,D-diff(f)]
    A\[exp(f[xm]),0.]    
end