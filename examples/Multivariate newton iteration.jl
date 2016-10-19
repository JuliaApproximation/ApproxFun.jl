using ApproxFun,DualNumbers,Plots
    import ApproxFun: DualFun, jacobian

x=Fun([-π,π])
N = (T,u) -> [u(-π)-1;u(π)-1;u'(π)-1;u'' - u - T*x/π - T^3*x^3/(π^3*6)]

T0=Tk=0.
    u0=uk=one(x)

    (T,u)=(T0,u0)
    for k=1:3
        J1=map(f->Fun(dualpart(coefficients(f)),space(f)),N(dual(T,1),u))
        J2=map(jacobian,N(T,DualFun(u)))
        J=[J1 J2]
        Tk,uk=J\N(T,u)
        T,u=T-Number(Tk),u-uk
    end


plot(u)
@show norm(N(T,u)[end]) # 1.39e-16
