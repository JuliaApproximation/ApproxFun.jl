
export Chebyshev


typealias Chebyshev Ultraspherical{0}


Space(d::Interval)=Chebyshev(d)
canonicalspace(S::Ultraspherical)=Chebyshev(domain(S))

function coefficients(g::Vector,::ConstantSpace,::Chebyshev)
    @assert length(g)==1
    g
end

function coefficients(g::Vector,::Chebyshev,::ConstantSpace)
    @assert length(g)==1
    g
end


## Transform

transform(::Chebyshev,vals::Vector,plan::Function)=chebyshevtransform(vals,plan)
itransform(::Chebyshev,cfs::Vector,plan::Function)=ichebyshevtransform(cfs,plan)
plan_transform(::Chebyshev,vals::Vector)=plan_chebyshevtransform(vals)
plan_itransform(::Chebyshev,cfs::Vector)=plan_ichebyshevtransform(cfs)

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


## Piecewise union

# union_rule dictates how to create a space that both spaces can be converted to
# in this case, it means
function union_rule{S<:Ultraspherical}(s1::PiecewiseSpace{S},s2::PiecewiseSpace{S})
    PiecewiseSpace(map(S,merge(domain(s1),domain(s2)).domains))
end

function union_rule{S<:Ultraspherical}(s1::PiecewiseSpace{S},s2::S)
    PiecewiseSpace(map(S,merge(domain(s1),domain(s2)).domains))
end



## Multivariate


#TODO: adaptive
for op in (:(Base.sin),:(Base.cos))
    @eval ($op){S<:Chebyshev,V<:Chebyshev}(f::ProductFun{S,V})=ProductFun(chebyshevtransform($op(values(f))),space(f))
end



reverseorientation(f::Fun{Chebyshev})=Fun(alternatesign!(copy(f.coefficients)),Chebyshev(reverse(domain(f))))
