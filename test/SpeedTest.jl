using ApproxFun, Base.Test


c = rand(1000)
x=rand(10000)
f=Fun(c,Chebyshev)
y=f(x)
y=f(x)

@time y=f(x)
println("Clenshaw large coeffs, many points: Time should be ~0.024")
# 0.012482274  with unsafe_view
# 0.024306262 with inbounds

y=f(.1)
y=f(.1)
y=f(.1)

@time y=f(.1);
println("Clenshaw large coeffs, 1 point: Time should be ~6e-6")

# @time is 8.853e-6 seconds


f=Fun(exp)
x=sample(f,100000)
x=sample(f,100000)
@time x=sample(f,100000)
println("Sample: Time should be ~0.25")
# 0.213793292 with unsafe_view
# 0.268162181 with inbounds


f=Fun(x->cos(x),20)
roots(f)
roots(f)
@time for k=1:100
    roots(f)
end
println("Small roots: Time should be ~0.015")



f=Fun(x->cos(1000x),1000)
roots(f)
roots(f)
@time roots(f)
println("Roots: Time should be ~0.13")
