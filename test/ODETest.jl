m=100.;
d=Interval(-m,5.);
Bm=EvaluationOperator(d.a,d);
Bp=EvaluationOperator(d.b,d);
B=[Bm,Bp];
D2=diff(d,2);
X=DifferentialOperator([Fun(x->x,d)],d);

u=[B,D2-X]\[airyai(d.a),airyai(d.b),0.];

@assert abs(u[0.]-airyai(0.)) < 10eps()
