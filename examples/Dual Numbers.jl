##
# Demo of using DualNumbers to differentiate with respect to
# changes in parameters
##




using DualNumbers

# We need to implement some functionality for the ApproxFun constructor to work
Base.real{T}(::Type{Dual{T}})=Dual{ApproxFun.real(T)}
Base.sinpi(x::Dual)=sin(π*x)

using ApproxFun
ApproxFun.plan_chebyshevtransform{D<:Dual}(v::Vector{D})=ApproxFun.plan_chebyshevtransform(real(v))
ApproxFun.chebyshevtransform{D<:Dual}(v::Vector{D},plan...)=dual(chebyshevtransform(real(v),plan...),chebyshevtransform(epsilon(v),plan...))
ApproxFun.chop!(f::Fun,d::Dual)=chop!(f,real(d))


# what is the derivative of the first coefficient with respect to the first argument of the domain?
f=Fun(exp,Interval(dual(1.0,1),dual(2.0)))
epsilon(f.coefficients[1])
h=0.00001;(Fun(exp,Interval(1.0+h,2.0)).coefficients[1]-Fun(exp,Interval(1.0,2.0)).coefficients[1])/h

