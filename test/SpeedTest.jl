using ApproxFun

c = rand(1000);
x=rand(10000);
f=Fun(c);

@time y=f[x]
#0.012482274


@ time y=f[.1]
#0.000352981


f=Fun(exp)

@time x=sample(f,100000)
#0.213793292


