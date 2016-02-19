using ApproxFun, Base.Test


## Dual Numbers


f=Fun(exp,Interval(dual(1.0,1),dual(2.0)),20)
@test_approx_eq Fun(h->Fun(exp,Interval(1.0+h,2.0)).coefficients[1],[0.,.1])'(0.) epsilon(f.coefficients[1])

