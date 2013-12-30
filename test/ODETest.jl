using ApproxFun


m=100.;
d=Interval(-m,5.);
Bm=EvaluationOperator(d.a,d);
Bp=EvaluationOperator(d.b,d);
B=[Bm,Bp];
D2=diff(d,2);
X=DifferentialOperator([Fun(x->x,d)],d);

u=[B,D2-X]\[airyai(d.a),airyai(d.b),0.];

@assert abs(u[0.]-airyai(0.)) < 10eps()



f=Fun(x->x.^2)
D=diff(f.domain)
@assert norm(D*f - diff(f)) < 10eps()


##Test versus exp

f=Fun(x->-x.^2)
g=Fun(t->exp(-t.^2))

@assert norm(Fun(t->exp(f[t]))-g)<10eps()

fp=diff(f);
Bm=EvaluationOperator(f.domain.a,f.domain);
u=[Bm,DifferentialOperator([-fp,1.],f.domain)]\[exp(f[f.domain.a]),0.];
@assert norm(u-g)<10eps()
