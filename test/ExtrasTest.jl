using ApproxFun, Base.Test


## Dual Numbers


f=Fun(exp,Interval(DualNumbers.dual(1.0,1),DualNumbers.dual(2.0)),20)
@test_approx_eq Fun(h->Fun(exp,Interval(1.0+h,2.0)).coefficients[1],[0.,.1])'(0.) DualNumbers.epsilon(f.coefficients[1])

