##
# Demo of using DualNumbers to differentiate with respect to
# changes in parameters
##


using DualNumbers, ApproxFun
# What is the derivative of the function (without differentiating the Chebyshev expansion)?
f=Fun(x->exp(dual(x,1)),-1..1)
dualpart(f)
# check versus differentiate
norm(realpart(f)'-dualpart(f))

# What is the derivative of the first coefficient with respect to the left endpoint?
f=Fun(exp,Segment(dual(1.0,1),dual(2.0)))
dualpart(f.coefficients[1])
# check versus finite difference calculation:
h=0.00001;(Fun(exp,Segment(1.0+h,2.0)).coefficients[1]-Fun(exp,Segment(1.0,2.0)).coefficients[1])/h

# Or an ApproxFun calculation:
Fun(h->Fun(exp,Segment(1.0+h,2.0)).coefficients[1],0..1)'(0.)

# What is the derivative of the first coefficient with respect to the exponential's constant?
f=Fun(x->exp(dual(x,x)),-1..1)
dualpart(f.coefficients[1])
# check versus finite difference calculation:
h=0.00001;(Fun(x->exp((1+h)x),-1..1).coefficients[1]-Fun(exp,-1..1).coefficients[1])/h
