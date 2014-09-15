include("bary.jl")
include("clenshaw.jl")
include("ultraspherical.jl")




##  Constructors



##TODO: No zero length funs
type Fun{T<:Union(Float64,Complex{Float64}),S<:FunctionSpace} 
    coefficients::Vector{T}
    space::S
end




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

##Coefficient routines
#TODO: domainscompatible?
coefficients(f::Fun,msp::FunctionSpace)=spaceconversion(f.coefficients,space(f),msp)
coefficients(f::Fun)=coefficients(f,ChebyshevSpace(domain(f)))

##Convert routines


Base.convert{T<:Number,S<:DomainSpace}(::Type{Fun{T,S}},x::Number)=x*ones(T,S(AnyDomain()))
Base.convert{T<:Number,S<:FunctionSpace}(::Type{Fun{Complex{Float64},S}},f::Fun{T,S})=Fun(convert(Vector{Complex{Float64}},f.coefficients),f.space)
Base.promote_rule{T<:Number,S<:FunctionSpace}(::Type{Fun{Complex{Float64},S}},::Type{Fun{T,S}})=Fun{Complex{Float64},S}
Base.promote_rule{T<:Number,IF<:Fun}(::Type{IF},::Type{T})=IF


Base.zero{T,S<:DomainSpace}(::Type{Fun{T,S}})=zeros(T,S(AnyDomain()))
Base.one{T,S<:DomainSpace}(::Type{Fun{T,S}})=ones(T,S(AnyDomain()))
for op in (:(Base.zeros),:(Base.ones))
    @eval ($op){T}(f::Fun{T})=$op(T,f.space)
end



## domain

for op = (:tocanonical,:tocanonicalD,:fromcanonical,:fromcanonicalD)
    @eval ($op)(f::Fun,x)=($op)(domain(f),x)
end

##Evaluation

evaluate{T,S}(f::Fun{T,S},x)=evaluate(Fun(f,domain(f)),x)  #Default is convert to Chebyshev
Base.getindex(f::Fun,x)=evaluate(f,x)



Base.first(f::Fun)=foldr(-,coefficients(f))
Base.last(f::Fun)=reduce(+,coefficients(f))


space(f::Fun)=f.space
spacescompatible(f::Fun,g::Fun)=spacescompatible(space(f),space(g))



##Data routines
values(f::Fun)=itransform(f.space,f.coefficients) 
points(f::Fun)=points(f.space,length(f))
Base.length(f::Fun)=length(f.coefficients)



## Manipulate length


pad!(f::Fun,n::Integer)=pad!(f.coefficients,n)
pad(f::Fun,n::Integer)=Fun(pad(f.coefficients,n),f.space)


function chop!{T}(f::Fun{T},tol::Real)
    chop!(f.coefficients,tol)
    if length(f.coefficients) == 0
        f.coefficients = [zero(T)]
    end
    
    f
end
chop(f::Fun,tol)=chop!(Fun(copy(f.coefficients),f.space),tol)
chop!(f::Fun)=chop!(f,eps())


## Addition and multiplication




for op = (:+,:-)
    @eval begin
        function ($op)(f::Fun,g::Fun)
            if spacescompatible(f,g)
                n = max(length(f),length(g))
                f2 = pad(f,n); g2 = pad(g,n)
            
                Fun(($op)(f2.coefficients,g2.coefficients),domain(f)!=AnyDomain()?f.space:g.space)
            else 
                $op(Fun(f,domain(f)),Fun(g,domain(g))) # convert to Chebyshev
            end
        end

        ($op){N<:Number,T<:Number}(f::Fun{T},c::N)=$op(f,c*ones(f))
        ($op){N<:Number,T<:Number}(c::N,f::Fun{T})=$op(c*ones(f),f)    
    end
end 



fasttimes(f2,g2)=Fun(chebyshevtransform(values(f2).*values(g2)),domain(f2))




for op = (:*,:.*,:./,:/)
    @eval ($op)(f::Fun,c::Number) = Fun(($op)(f.coefficients,c),f.space)
end 

-(f::Fun)=Fun(-f.coefficients,f.space)
-(c::Number,f::Fun)=-(f-c)


for op = (:*,:.*,:+)
    @eval ($op)(c::Number,f::Fun)=($op)(f,c)
end




function .^(f::Fun,k::Integer)
    if k == 0
        1.
    elseif k > 0
        f.*f.^(k-1)
    else
        f./f.^(k+1)
    end
end


# multiplying 2 Funs, we assume this can be done by transform
# the parametrizations are to make it the broadest definition

function .*{T,N,S}(f::Fun{T,S},g::Fun{N,S})
    @assert domainscompatible(f,g)
    #TODO Coefficient space version
    n = length(f) + length(g) - 1
    f2 = pad(f,n); g2 = pad(g,n)
    
    sp=space(f)
    chop!(Fun(transform(sp,values(f2).*values(g2)),sp),10eps())
end



## Norm

import Base.norm

norm(f::Fun)=real(sqrt(sum(f.*conj(f))))



## Mapped functions

import Base.imag, Base.real, Base.conj

for op = (:real,:imag,:conj) 
    ##TODO: this assumes real space
    @eval ($op)(f::Fun) = Fun(($op)(f.coefficients),f.space)
end

Base.abs2(f::Fun{Float64})=f.^2
Base.abs2(f::Fun{Complex{Float64}})=real(f).^2+imag(f).^2

##  integration

function Base.cumsum(f::Fun)
    cf = integrate(f)
    cf - first(cf)
end

function Base.sum{T}(f::Fun{T})
    cf=integrate(f)
    last(cf) - first(cf)
end




function differentiate(f::Fun,k::Integer)
    @assert k >= 0
    (k==0)?f:differentiate(differentiate(f),k-1)
end

Base.diff(f::Fun,n...)=differentiate(f,n...)



==(f::Fun,g::Fun) =  (f.coefficients == g.coefficients && f.space == g.space)



