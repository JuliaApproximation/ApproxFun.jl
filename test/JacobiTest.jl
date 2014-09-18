using ApproxFun, Base.Test

f=Fun(exp,JacobiSpace(2.,.5))
@test_approx_eq f[.1] exp(.1)

f=Fun(x->cos(100x),JacobiSpace(2.124,.5),500)
@test_approx_eq f[.1] cos(100*.1)
