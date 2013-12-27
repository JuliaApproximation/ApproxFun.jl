

function Base.exp(f::IFun)
    fp=diff(f);
    Bm=EvaluationOperator(f.domain.a,f.domain);
    u=[Bm,DifferentialOperator([-fp,1.],f.domain)]\[exp(f[f.domain.a]),0.];
end