using ApproxFun
setplotter("PyPlot")

#undamped simple harmonic oscillator
#f(x,y) = [ y(2), -y(1) ]

#undamped pendulum
f(x,y) = [ y[2], -sin(y[1]) ]

#van der Pol oscillator. reaches coefficients length limit after a few iterations.
#μ = 1.
#f(x,y) = [ y[2], μ*(1-y[1]^2)*y[2] - y[1] ]


y0 = [3., 0.]

D = Interval(0.,20.)

xt = Fun(identity, D)

y0 = ApproxFun.devec(Fun[Fun(x->y, D, method="abszerocoefficients") for y in y0])
yn = y0
ApproxFun.plot(vec(yn)[1])


while true
	err = sum([maximum(fun^2) for fun in vec( diff(yn) -  ApproxFun.devec(f(xt,vec(yn))) )])
	println("Length: $(length(yn)) \tMaximum point-wise squared error: $err")
	if err < 1e-10
		break
	end

	ApproxFun.plot(vec(yn)[1])
	integrand = Fun(x::Float64->f(x,yn(x)), D)
	yn = cumsum(integrand) + y0
end

xts = D.a:0.01:D.b
PyPlot.plot(xts,vec(yn)[1]([xts]);lw=3)
PyPlot.plot(xts,vec(yn)[2]([xts]);lw=3)
