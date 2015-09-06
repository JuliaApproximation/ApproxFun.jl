export newton

immutable DualFun{F,T}
    f::F
    J::T
end
DualFun(f::Fun)=DualFun(f,
                        SpaceOperator(IdentityOperator(),
                              space(f),
                              space(f)))

differentiate(d::DualFun)=DualFun(d.f',Derivative()*d.J)
Base.transpose(d::DualFun)=differentiate(d)

^(d::DualFun,k::Integer)=DualFun(d.f^k,k*d.f^(k-1)*d.J)
+(a::DualFun,b::DualFun)=DualFun(a.f+b.f,a.J+b.J)
for OP in (:+,:-)
    @eval $OP(a::DualFun,b::Number)=DualFun($OP(a.f,b),a.J)
end
*(f::Number,d::DualFun)=DualFun(f*d.f,f*d.J)
*(f::Fun,d::DualFun)=DualFun(f*d.f,f*d.J)
*(a::DualFun,b::DualFun)=DualFun(a.f*b.f,a.f*b.J+b.f*a.J)

Base.getindex(d::DualFun,x)=DualFun(d.f[x],Evaluation(rangespace(d.J),x)*d.J)

jacobian(d::DualFun)=d.J

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
