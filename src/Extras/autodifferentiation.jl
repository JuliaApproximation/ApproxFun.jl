export newton

immutable DualFun{F,T}
    f::F
    J::T
end
DualFun(f::Fun)=DualFun(f,
                        SpaceOperator(IdentityOperator(),
                              space(f),
                              space(f)))


domain(df::DualFun)=domain(df.f)

differentiate(d::DualFun)=DualFun(d.f',Derivative(rangespace(d.J))*d.J)
integrate(d::DualFun)=DualFun(integrate(d.f),Integral(rangespace(d.J))*d.J)
function Base.cumsum(d::DualFun)
    Q=Integral(rangespace(d.J))*d.J
    DualFun(cumsum(d.f),(I-Evaluation(rangespace(Q),false))*Q)
end


Base.transpose(d::DualFun)=differentiate(d)

^(d::DualFun,k::Integer)=DualFun(d.f^k,k*d.f^(k-1)*d.J)

for OP in (:+,:-)
    @eval begin
        $OP(a::DualFun,b::Number)=DualFun($OP(a.f,b),a.J)
        $OP(a::DualFun,b::DualFun)=DualFun($OP(a.f,b.f),$OP(a.J,b.J))
    end
end
-(a::DualFun)=DualFun(-a.f,-a.J)

*(f::Number,d::DualFun)=DualFun(f*d.f,f*d.J)
*(f::Fun,d::DualFun)=DualFun(f*d.f,f*d.J)
*(a::DualFun,b::DualFun)=DualFun(a.f*b.f,a.f*b.J+b.f*a.J)

Base.call(d::DualFun,x)=DualFun(d.f(x),Evaluation(rangespace(d.J),x)*d.J)
Base.first(d::DualFun)=DualFun(first(d.f),Evaluation(rangespace(d.J),false)*d.J)
Base.last(d::DualFun)=DualFun(last(d.f),Evaluation(rangespace(d.J),true)*d.J)

jacobian(d::DualFun)=d.J




function Operator(f,ds::Space)
    if (isgeneric(f)&&applicable(f,0)) || (!isgeneric(f)&&arglength(f)==1)
        df=f(DualFun(Fun(ds)))
    elseif (isgeneric(f)&&applicable(f,0,0)) || (!isgeneric(f)&&arglength(f)==2)
        df=f(Fun(ds),DualFun(Fun(ds)))
    else
        error("Not implemented")
    end

    if isa(df,Array)
        map(u->u.J,df)
    else
        df.J
    end
end

Operator(f,d)=Operator(f,Space(d))
Operator(f)=Operator(f,Chebyshev())  #TODO: UnsetSpace



# full operator should be
# N=u->[B*u-bcs;...]
function newton(N,u0;maxiterations=100,tolerance=1E-15)
    u=u0
    for k=1:maxiterations
        DF=N(DualFun(u))
        J=map(jacobian,DF)
        F=map(d->d.f,DF)
        unew=u-J\F
        if norm(unew-u)≤10tolerance
            return unew
        else
            u=chop(unew,tolerance)
        end
    end
end
