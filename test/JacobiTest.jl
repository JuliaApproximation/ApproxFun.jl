using ApproxFun, Base.Test

f=Fun(exp,JacobiSpace(2.,.5))
@test_approx_eq f[.1] exp(.1)

f=Fun(x->cos(100x),JacobiSpace(2.124,.5),500)
@test_approx_eq f[.1] cos(100*.1)


sp=JacobiSpace(2.124,.5)
f=Fun(exp,sp)
sp2=JacobiSpace(2.124,1.5)
f2=Fun(exp,sp2)
sp3=JacobiSpace(3.124,1.5)
f3=Fun(exp,sp3)
sp4=JacobiSpace(4.124,2.5)
f4=Fun(exp,sp4)
@test norm((Fun(f,sp2)-f2).coefficients)<10eps()
@test norm((Fun(f,sp3)-f3).coefficients)<10eps()
@test norm((Fun(f,sp4)-f4).coefficients)<10eps()