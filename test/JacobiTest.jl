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


m=20
f=Fun(x->((1-x)/2).^m.*exp(x),JacobiWeightSpace(0.,m,JacobiSpace(2m+1,0.)))
@test abs(f[.1]-(x->((1-x)/2).^m.*exp(x))(.1))<10eps()


m=10
f=Fun(x->besselj(m,m*(1-x)),JacobiWeightSpace(0.,m,JacobiSpace(2m+1,0.)))
@test_approx_eq f[0.] besselj(m,m)



## Disk

f=(x,y)->exp(x.*sin(y))
u=ProductFun(f,Disk(),50,51)
@test_approx_eq u[.1,.1] f(.1,.1)

