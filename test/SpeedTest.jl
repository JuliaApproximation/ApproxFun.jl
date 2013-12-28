using ApproxFun

c = rand(1000);
x=rand(10000);
f=Fun(c);

tic();y=f[x];toc()
#0.012482274


tic();y=f[.1];toc()
#0.000352981


f=Fun(exp)

@time x=sample(f,10000)
#0.015744034 or 0.039095687


