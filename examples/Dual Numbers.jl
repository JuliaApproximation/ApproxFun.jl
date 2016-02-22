##
# Demo of using DualNumbers to differentiate with respect to
# changes in parameters
##


using DualNumbers, ApproxFun

# what is the derivative of the first coefficient with respect to the first argument of the domain?
f=Fun(exp,Interval(dual(1.0,1),dual(2.0)),20)
epsilon(f.coefficients[1])

# check versus finite difference calculation:
h=0.00001;(Fun(exp,Interval(1.0+h,2.0)).coefficients[1]-Fun(exp,Interval(1.0,2.0)).coefficients[1])/h

# Or an ApproxFun calculation:
Fun(h->Fun(exp,Interval(1.0+h,2.0)).coefficients[1],[0.,.1])'(0.)
