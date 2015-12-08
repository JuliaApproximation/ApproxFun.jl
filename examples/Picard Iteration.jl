using Plots,ApproxFun

#undamped simple harmonic oscillator
#f =(x,y)-> [ y(2), -y(1) ]

#undamped pendulum
f =(x,y) -> [ y[2], -sin(y[1]) ]

#van der Pol oscillator. reaches coefficients length limit after a few iterations.
#μ = 1.
#f(x,y) = [ y[2], μ*(1-y[1]^2)*y[2] - y[1] ]


y0 = [3., 0.]
D = Interval(0.,20.)
xt = Fun(identity, D)
y0=Fun(x->y0,D)
yn = y0
plot(yn[1])

while true
	err = sum([maximum(fun^2) for fun in vec(yn'-f(x,yn))])
	println("Length: $(length(yn)) \tMaximum point-wise squared error: $err")
	if err < 1e-10
		break
	end

    plot(yn[1])
	integrand = Fun(f(x,yn))
    yn = chop(cumsum(integrand) + y0,100eps())
end


xts = D.a:0.01:D.b
plot(xts,yn[1])
plot!(xts,yn[2])
