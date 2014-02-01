using ApproxFun

c = rand(1000);
x=rand(10000);
f=Fun(c);

@time y=f[x]
# 0.012482274  with unsafe_view
# 0.024306262 with inbounds


@ time y=f[.1]
# 8.853e-6 seconds


f=Fun(exp)

@time x=sample(f,100000)
# 0.213793292 with unsafe_view
# 0.268162181 with inbounds


