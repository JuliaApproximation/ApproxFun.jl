
export Chebyshev


typealias Chebyshev Ultraspherical{0}


Space(d::IntervalDomain)=Chebyshev(d)
canonicalspace(S::Ultraspherical)=Chebyshev(domain(S))

function spaceconversion(g::Vector,::ConstantSpace,::Chebyshev)
    @assert length(g)==1
    g
end

function spaceconversion(g::Vector,::Chebyshev,::ConstantSpace)
    @assert length(g)==1
    g
end


## Transform

transform(::Chebyshev,vals::Vector)=chebyshevtransform(vals)
itransform(::Chebyshev,cfs::Vector)=ichebyshevtransform(cfs)


## Evaluation

evaluate(f::Fun{Chebyshev},x)=clenshaw(f.coefficients,tocanonical(f,x))

## Calculus


# diff T -> U, then convert U -> T
integrate(f::Fun{Chebyshev})=Fun(chebyshevintegrate(domain(f),f.coefficients),f.space)
chebyshevintegrate(d::Interval,cfs::Vector)=fromcanonicalD(d,0)*ultraint!(ultraconversion(cfs))   


differentiate(f::Fun{Chebyshev})=Fun(chebyshevdifferentiate(domain(f),f.coefficients),f.space)
chebyshevdifferentiate(d::Interval,cfs::Vector)=tocanonicalD(d,0)*ultraiconversion(ultradiff(cfs))
chebyshevdifferentiate(d::IntervalDomain,cfs::Vector)=(Fun(x->tocanonicalD(d,x),d).*Fun(differentiate(Fun(cfs)),d)).coefficients


## identity_fun

identity_fun(d::Chebyshev)=identity_fun(domain(d))



## 2D fast values

function ApproxFun.values{T}(f::TensorFun{Chebyshev,Chebyshev,T})
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



## Piecewise union

# union_rule dictates how to create a space that both spaces can be converted to
# in this case, it means 
function union_rule{S<:Ultraspherical}(s1::PiecewiseSpace{S},s2::PiecewiseSpace{S})
    PiecewiseSpace(map(S,merge(domain(s1),domain(s2)).domains))
end

function union_rule{S<:Ultraspherical}(s1::PiecewiseSpace{S},s2::S)
    PiecewiseSpace(map(S,merge(domain(s1),domain(s2)).domains))
end



## Multivaraite


#TODO: adaptive
for op in (:(Base.sin),:(Base.cos))
    @eval ($op){S<:Chebyshev,V<:Chebyshev,T}(f::TensorFun{S,V,T})=TensorFun(chebyshevtransform($op(values(f))),space(f))
end

