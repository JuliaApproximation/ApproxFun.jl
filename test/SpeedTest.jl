using ApproxFun, Base.Test

c = rand(1000)
x=rand(10000)
f=Fun(c)
y=f[x]
y=f[x]

@time y=f[x]
println("Time should be ~0.024")
# 0.012482274  with unsafe_view
# 0.024306262 with inbounds

y=f[.1]
y=f[.1]
y=f[.1]

@time y=f[.1];
println("Time should be ~9e-6")

# @time is 8.853e-6 seconds


f=Fun(exp)
x=sample(f,100000)
x=sample(f,100000)
@time x=sample(f,100000)
println("Time should be ~0.27")
# 0.213793292 with unsafe_view
# 0.268162181 with inbounds


