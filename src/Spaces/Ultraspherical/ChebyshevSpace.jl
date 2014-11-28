

typealias ChebyshevSpace UltrasphericalSpace{0}


Space(d::IntervalDomain)=ChebyshevSpace(d)
canonicalspace(S::UltrasphericalSpace)=ChebyshevSpace(domain(S))

function spaceconversion(g::Vector,::ConstantSpace,::ChebyshevSpace)
    @assert length(g)==1
    g
end

function spaceconversion(g::Vector,::ChebyshevSpace,::ConstantSpace)
    @assert length(g)==1
    g
end


## Transform

transform(::ChebyshevSpace,vals::Vector)=chebyshevtransform(vals)
itransform(::ChebyshevSpace,cfs::Vector)=ichebyshevtransform(cfs)


## Evaluation

evaluate(f::Fun{ChebyshevSpace},x)=clenshaw(f.coefficients,tocanonical(f,x))

## Calculus


# diff T -> U, then convert U -> T
integrate(f::Fun{ChebyshevSpace})=Fun(chebyshevintegrate(domain(f),f.coefficients),f.space)
chebyshevintegrate(d::Interval,cfs::Vector)=fromcanonicalD(d,0)*ultraint(ultraconversion(cfs))   


differentiate(f::Fun{ChebyshevSpace})=Fun(chebyshevdifferentiate(domain(f),f.coefficients),f.space)
chebyshevdifferentiate(d::Interval,cfs::Vector)=tocanonicalD(d,0)*ultraiconversion(ultradiff(cfs))
chebyshevdifferentiate(d::IntervalDomain,cfs::Vector)=(Fun(x->tocanonicalD(d,x),d).*Fun(diff(Fun(cfs)),d)).coefficients


## identity_fun

identity_fun(d::ChebyshevSpace)=identity_fun(domain(d))



## 2D fast values

function ApproxFun.values{T}(f::TensorFun{ChebyshevSpace,ChebyshevSpace,T})
    n,m=size(f)
    M=Array(T,n,m)
    f1=pad(f.coefficients[1].coefficients,n)
    planc=plan_chebyshevtransform(f1)
    M[:,1]=ichebyshevtransform(f1,planc)
    for k=2:m
        M[:,k]=ichebyshevtransform(pad(f.coefficients[k].coefficients,n),planc)
    end
    f2=vec(M[1,:])
    planr=plan_chebyshevtransform(f2)
    M[1,:]=ichebyshevtransform(f2,planr)
    for k=2:n
        M[k,:]=ichebyshevtransform(vec(M[k,:]),planr)
    end
    
    M
end