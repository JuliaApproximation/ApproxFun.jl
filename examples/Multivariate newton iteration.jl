using ApproxFun,DualNumbers,Plots
    import ApproxFun: DualFun, jacobian

x=Fun(-π..π)
N = (T,u) -> [u(-π)-1;u(π)-1;u'(π)-1;u'' - u - sin(T*x/π)]

T0=Tk=0.
    u0=uk=one(x)

    (T,u)=(T0,u0)
    for k=1:4
        J1=jacobian.(N(DualFun(Fun(T)),u))
        J2=jacobian.(N(T,DualFun(u)))
        J=Operator[J1 J2]
        Tk,uk=J\N(T,u)
        @show Tk
        T,u=T-Number(Tk),u-uk
    end


plot(u)
norm(N(T,u)[end])  #9e-17
