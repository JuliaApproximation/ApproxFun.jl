##
# Demo of using DualNumbers to differentiate with respect to
# changes in parameters
##

using ApproxFun, DualNumbers

# What is the derivative of the function (without differentiating the Chebyshev expansion)?
f=Fun(x->exp(Dual(x,1)),[-1,1])
epsilon(f)
# check versus differentiate
norm(value(f)'-epsilon(f))

# What is the derivative of the first coefficient with respect to the left endpoint?
f=Fun(exp,Interval(Dual(1.0,1),Dual(2.0)))
epsilon(f.coefficients[1])
# check versus finite difference calculation:
h=0.00001;(Fun(exp,Interval(1.0+h,2.0)).coefficients[1]-Fun(exp,Interval(1.0,2.0)).coefficients[1])/h

# What is the derivative of the first coefficient with respect to the exponential's constant?
f=Fun(x->exp(Dual(x,x)),[-1,1])
epsilon(f.coefficients[1])
# check versus finite difference calculation:
h=0.00001;(Fun(x->exp((1+h)x),[-1,1]).coefficients[1]-Fun(exp,[-1,1]).coefficients[1])/h

