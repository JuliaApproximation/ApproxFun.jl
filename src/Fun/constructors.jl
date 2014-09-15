Fun{T<:Union(Int64,Complex{Int64})}(coefs::Vector{T},d::FunctionSpace)=Fun(1.0*coefs,d)


function Fun(f::Function,d::DomainSpace,n::Integer)
    pts=points(d,n)
    f1=f(pts[1])
    T=typeof(f1)
        
    if T <: Vector
        Fun{typeof(f1[1]),typeof(d)}[Fun(x->f(x)[k],d,n) for k=1:length(f1)]
    elseif T <: Array
        Fun{typeof(f1[1,1]),typeof(d)}[Fun(x->f(x)[k,j],d,n) for k=1:size(f1,1),j=1:size(f1,2)]    
    else
        vals=T[f(x) for x in pts]
        Fun(transform(d,vals),d)
    end
end

# the following is to avoid ambiguity
for SS in (:IntervalDomainSpace,:PeriodicDomainSpace)
    @eval Fun(f::Fun,d::($SS))=Fun(coefficients(f,d),d)
    @eval Fun{T<:($SS)}(f::Fun,::Type{T})=Fun(f,T(domain(f)))
    @eval Fun{T<:($SS)}(c::Number,::Type{T})=Fun(c,T(AnyDomain()))
end


Fun{T<:IntervalDomainSpace}(f,::Type{T})=Fun(f,T(Interval()))
Fun{T<:PeriodicDomainSpace}(f,::Type{T})=Fun(f,T(PeriodicInterval()))
Fun{T<:IntervalDomainSpace}(f,::Type{T},n::Integer)=Fun(f,T(Interval()),n)
Fun{T<:PeriodicDomainSpace}(f,::Type{T},n::Integer)=Fun(f,T(PeriodicInterval()),n)

Fun(f,d::IntervalDomain)=Fun(f,ChebyshevSpace(d))
Fun(f,d::IntervalDomain,n)=Fun(f,ChebyshevSpace(d),n)
Fun(f,d::PeriodicDomain)=Fun(f,LaurentSpace(d))
Fun(f,d::PeriodicDomain,n)=Fun(f,LaurentSpace(d),n)

Fun(f::Function,n::Integer)=Fun(f,Interval(),n)
Fun{T<:Number}(f::Function,d::Vector{T},n::Integer)=Fun(f,Interval(d),n)
Fun(cfs::Vector)=Fun(1.0*cfs,Interval())
Fun{T<:Number}(cfs::Vector,d::Vector{T})=Fun(1.0*cfs,Interval(d))
Fun(f::Function)=Fun(f,Interval())
Fun{T<:Number}(f::Function,d::Vector{T})=Fun(f,Interval(d))


Fun{T<:Number}(f::Fun,d::Vector{T})=Fun(coefficients(f),d)
Fun(f::Fun)=Fun(coefficients(f))

Fun(c::Number)=Fun([c])

for DD in (:IntervalDomain,:PeriodicDomain)
    @eval Fun(c::Number,d::($DD))=Fun([c],d)
end

Fun(c::Number,d)=Fun([c],d)

## List constructor

Fun{T<:Domain}(c::Number,dl::Vector{T})=map(d->Fun(c,d),dl)
Fun{T<:Domain}(f,dl::Vector{T})=map(d->Fun(f,d),dl)

## Adaptive constructors

function randomFun(f::Function,d::IntervalDomain)
    @assert d == Interval()

    #TODO: implement other domains
    
    Fun(chebyshevtransform(randomadaptivebary(f)),d)
end


function veczerocfsFun(f::Function,d::IntervalDomain)
    #reuse function values

    tol = 200*eps()

    for logn = 4:20
        cf = Fun(f, d, 2^logn + 1)
        cfs=coefficients(cf)  ##TODO: general domain
        
        if norm(cfs[:,end-8:end],Inf) < tol*norm(cfs[:,1:8],Inf)
            nrm=norm(cfs,Inf)
            return map!(g->chop!(g,10eps()*nrm),cf)
        end
    end
    
    warn("Maximum length reached")
    
    Fun(f,d,2^21 + 1)
end

function zerocfsFun(f::Function,d::DomainSpace)
    #reuse function values

    if isa(f(fromcanonical(d,0.)),Vector)       #TODO: is 0. going to always be in canonical?
        return veczerocfsFun(f,domain(d))
    end

    tol = 200*eps()

    for logn = 4:20
        cf = Fun(f, d, 2^logn + 1)
        
        if maximum(abs(cf.coefficients[end-8:end])) < tol*maximum(abs(cf.coefficients[1:8]))
            return chop!(cf,10eps()*maximum(abs(cf.coefficients)))
        end
    end
    
    warn("Maximum length reached")
    
    Fun(f,d,2^21 + 1)
end




function abszerocfsFun(f::Function,d::DomainSpace)
    #reuse function values

    tol = 200eps();

    for logn = 4:20
        cf = Fun(f, d, 2^logn + 1)
        
        if maximum(abs(cf.coefficients[end-8:end])) < tol
            return chop!(cf,10eps())
        end
    end
    
    warn("Maximum length reached")
    
    Fun(f,d,2^21 + 1)
end


function Fun(f::Function, d::DomainSpace; method="zerocoefficients")
    if f==identity
        identity_fun(d)
    elseif f==zero
        zeros(Float64,d)
    elseif f==one
        ones(Float64,d)
    elseif method == "zerocoefficients"
        zerocfsFun(f,d)
    elseif method == "abszerocoefficients"
        abszerocfsFun(f,d)
    else
        randomFun(f,d)    
    end
end
Fun(f::Function,d::IntervalDomain;opts...)=Fun(f,ChebyshevSpace(d);opts...)





## Aliases


FFun(x,d::PeriodicDomain)=Fun(x,LaurentSpace(d))
FFun(x,d::PeriodicDomain,n...)=Fun(x,LaurentSpace(d),n...)
FFun{T<:Number}(x,d::Vector{T})=Fun(x,LaurentSpace(d))
FFun{T<:Number}(x,d::Vector{T},n...)=Fun(x,LaurentSpace(d),n...)
FFun(f,n::Integer)=Fun(f,LaurentSpace(PeriodicInterval()),n)
FFun(f)=Fun(f,LaurentSpace(PeriodicInterval()))
