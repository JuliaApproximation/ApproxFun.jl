using ApproxFun, Base.Test
using DualNumbers

f = Fun(sin,Interval(-6.,6.))

println("Automatic differentiation and ApproxFun Differentiation comparison")

#this is a good test to ensure that Funs can be evaluated on non-real, non-complex number types like "Dual".
df_auto = Fun(x->epsilon( f[Dual(x,1.0)] ), Interval(-6.,6.))
df_poly = diff(f)

err = df_auto - df_poly

@test maximum(abs(err.coefficients)) < 10eps()	

