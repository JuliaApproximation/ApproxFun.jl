include("bary.jl")
include("clenshaw.jl")
include("ultraspherical.jl")

import Base.chop

export plan_chebyshevtransform

## transforms


function negateeven!(x)
    for k =2:2:length(x)
        x[k] *= -1.
    end
    
    x
end

plan_chebyshevtransform(x)=length(x)==1?identity:FFTW.plan_r2r(x, FFTW.REDFT00)
chebyshevtransform(x)=chebyshevtransform(x,plan_chebyshevtransform(x))
function chebyshevtransform(x::Vector,plan::Function)
    if(length(x) == 1)
        x
    else
        ret = plan(x)
        ret[1] *= .5
        ret[end] *= .5   
        negateeven!(ret)
        ret*=1./(length(ret)-1)
        
        ret
    end
end


ichebyshevtransform(x)=ichebyshevtransform(x,plan_chebyshevtransform(x))
function ichebyshevtransform(x::Vector,plan::Function)
    if(length(x) == 1)
        x
    else
        ##TODO: make thread safe
        x[1] *= 2.;
        x[end] *= 2.;
        
        ret = chebyshevtransform(x,plan)::typeof(x)
        
        x[1] *= .5;
        x[end] *= .5;
        
        ret[1] *= 2.;
        ret[end] *= 2.;
        
        negateeven!(ret)
        
        ret *= .5*(length(x) - 1)
        
        flipud(ret)
    end
end





##  Constructors




type IFun{T<:Union(Float64,Complex{Float64}),D<:IntervalDomain} <: AbstractFun
    coefficients::Vector{T}
    domain::D
end

function IFun(f::Function,d::Domain,n::Integer)
    pts=points(d,n)
    f1=f(pts[1])
    T=typeof(f1)
        
    if T <: Vector
        IFun{typeof(f1[1]),typeof(d)}[IFun(x->f(x)[k],d,n) for k=1:length(f1)]
    elseif T <: Array
        IFun{typeof(f1[1,1]),typeof(d)}[IFun(x->f(x)[k,j],d,n) for k=1:size(f1,1),j=1:size(f1,2)]    
    else
        vals=T[f(x) for x in pts]
        IFun(chebyshevtransform(vals),d)
    end
end


IFun(f::Function,n::Integer)=IFun(f,Interval(),n)
IFun{T<:Number}(f::Function,d::Vector{T},n::Integer)=IFun(f,Interval(d),n)
IFun(cfs::Vector)=IFun(1.0*cfs,Interval())
IFun{T<:Number}(cfs::Vector,d::Vector{T})=IFun(1.0*cfs,Interval(d))
IFun(cfs::Vector,d::IntervalDomain)=IFun(1.0*cfs,d)
IFun(f::Function)=IFun(f,Interval())
IFun{T<:Number}(f::Function,d::Vector{T})=IFun(f,Interval(d))

IFun(f::IFun,d::IntervalDomain)=IFun(f.coefficients,d)
IFun{T<:Number}(f::IFun,d::Vector{T})=IFun(f.coefficients,d)
IFun(f::IFun)=IFun(f.coefficients)

IFun(c::Number)=IFun([c])
IFun(c::Number,d)=IFun([c],d)

## List constructor

IFun{T<:IntervalDomain}(c::Number,dl::Vector{T})=map(d->IFun(c,d),dl)
IFun{T<:IntervalDomain}(f,dl::Vector{T})=map(d->IFun(f,d),dl)

## Adaptive constructors

function randomIFun(f::Function,d::Domain)
    @assert d == Interval()

    #TODO: implement other domains
    
    IFun(chebyshevtransform(randomadaptivebary(f)),d)
end


function veczerocfsIFun(f::Function,d::Domain)
    #reuse function values

    tol = 200*eps()

    for logn = 4:20
        cf = IFun(f, d, 2^logn + 1)
        cfs=coefficients(cf)
        
        if norm(cfs[:,end-8:end],Inf) < tol*norm(cfs[:,1:8],Inf)
            nrm=norm(cfs,Inf)
            return map!(g->chop!(g,10eps()*nrm),cf)
        end
    end
    
    warn("Maximum length reached")
    
    IFun(f,d,2^21 + 1)
end

function zerocfsIFun(f::Function,d::Domain)
    #reuse function values

    if isa(f(first(d)),Vector)
        return veczerocfsIFun(f,d)
    end

    tol = 200*eps()

    for logn = 4:20
        cf = IFun(f, d, 2^logn + 1)
        
        if maximum(abs(cf.coefficients[end-8:end])) < tol*maximum(abs(cf.coefficients[1:8]))
            return chop!(cf,10eps()*maximum(abs(cf.coefficients)))
        end
    end
    
    warn("Maximum length reached")
    
    IFun(f,d,2^21 + 1)
end




function abszerocfsIFun(f::Function,d::Domain)
    #reuse function values

    tol = 200eps();

    for logn = 4:20
        cf = IFun(f, d, 2^logn + 1);
        
        if maximum(abs(cf.coefficients[end-8:end])) < tol
            return chop!(cf,10eps())
        end
    end
    
    warn("Maximum length reached")
    
    IFun(f,d,2^21 + 1)
end


function IFun(f::Function, d::Domain; method="zerocoefficients")
    if f==identity
        identity_fun(d)
    elseif f==zero
        IFun([0.0],d)
    elseif f==one
        IFun([1.0],d)    
    elseif method == "zerocoefficients"
        zerocfsIFun(f,d)
    elseif method == "abszerocoefficients"
        abszerocfsIFun(f,d)
    else
        randomIFun(f,d)    
    end
end

##Coefficient routines

coefficients(f::IFun)=f.coefficients
coefficients(f::IFun,m::Integer)=ultraconversion(f.coefficients,m)

##Convert routines
Base.convert{T<:Number,D<:IntervalDomain}(::Type{IFun{T,D}},x::Number)=IFun([1.*x])
Base.convert(::Type{IFun},x::Int64)=IFun([1.*x])
Base.convert{D<:IntervalDomain}(::Type{IFun{Float64,D}},x::IFun)=1.*x

##Evaluation


Base.getindex(f::IFun,x)=evaluate(f,x)
evaluate(f::IFun,x)=clenshaw(f.coefficients,tocanonical(f.domain,x))


Base.first(f::IFun)=foldr(-,f.coefficients)
Base.last(f::IFun)=reduce(+,f.coefficients)




##Data routines
values(f::IFun)=ichebyshevtransform(f.coefficients) 

points(f::IFun)=points(f.domain,length(f))


Base.length(f::IFun)=length(f.coefficients)



## Manipulate length


pad!(f::IFun,n::Integer)=pad!(f.coefficients,n)
pad(f::IFun,n::Integer)=IFun(pad(f.coefficients,n),f.domain)


function chop!(f::IFun,tol::Real)
    chop!(f.coefficients,tol)
    if length(f.coefficients) == 0
        f.coefficients = [0.]
    end
    
    f
end
chop(f::Union(IFun,Vector),tol)=chop!(deepcopy(f),tol)

chop!(f::Union(IFun,Vector))=chop!(f,eps())


## Addition and multiplication




for op = (:+,:-)
    @eval begin
        function ($op)(f::IFun,g::IFun)
            @assert f.domain == g.domain
        
            n = max(length(f),length(g))
            f2 = pad(f,n);
            g2 = pad(g,n);
            
            IFun(($op)(f2.coefficients,g2.coefficients),f.domain)
        end

        function ($op)(f::IFun,c::Number)
            f2 = deepcopy(f);
            
            f2.coefficients[1] = ($op)(f2.coefficients[1],c);
            
            f2
        end
    end
end 



function .*(f::IFun,g::IFun)
    @assert f.domain == g.domain
    #TODO Coefficient space version
    n = length(f) + length(g) - 1;
    f2 = pad(f,n);
    g2 = pad(g,n);
    
    chop!(IFun(chebyshevtransform(values(f2).*values(g2)),f.domain),10eps())
end

fasttimes(f2,g2)=IFun(chebyshevtransform(values(f2).*values(g2)),f2.domain)




for op = (:*,:.*,:./,:/)
    @eval ($op)(f::IFun,c::Number) = IFun(($op)(f.coefficients,c),f.domain)
end 

-(f::IFun)=IFun(-f.coefficients,f.domain)
-(c::Number,f::IFun)=-(f-c)


for op = (:*,:.*,:+)
    @eval ($op)(c::Number,f::IFun)=($op)(f,c)
end




function .^(f::IFun,k::Integer)
    if k == 0
        1.
    elseif k > 0
        f.*f.^(k-1)
    else
        f./f.^(k+1)
    end
end


## Norm

import Base.norm

norm(f::IFun)=real(sqrt(sum(f.*conj(f))))



## Mapped functions

import Base.imag, Base.real, Base.conj

for op = (:real,:imag,:conj) 
    @eval ($op)(f::IFun) = IFun(($op)(f.coefficients),f.domain)
end


## Differentiation and integration



function Base.cumsum(f::IFun)
    cf = integrate(f)
    cf - first(cf)
end

function Base.sum(f::IFun)
    cf=integrate(f)
    last(cf) - first(cf)
end



==(f::IFun,g::IFun) =  (f.coefficients == g.coefficients && f.domain == g.domain)



## Root finding


function complexroots(cin::Vector)
    c=chop(cin,10eps())
    if c == [] || length(c) == 1
        return []
    elseif length(c) == 2
        return [-c[1]/c[2]]
    else 
        n=length(c)-1;
        
        I = [ones(Int64,n),2:n-1,2:n];
        J=[1:n,3:n,1:n-1];
        V = [-c[end-1]/(2c[end]),.5-c[end-2]/(2c[end]),-c[end-3:-1:1]/(2c[end]),.5*ones(n-2),.5*ones(n-2),1];
        A=full(sparse(I,J,V));
        
        return eigvals(A)
    end
end


function complexroots(f::IFun)
    fromcanonical(f,complexroots(f.coefficients))
end



function roots(f::IFun)
    irts=map(real,filter!(x->abs(x)<=1.+10eps(),filter!(isreal,complexroots(f.coefficients))))
    
    map!(x->x>1.?1.:x,irts)
    map!(x->x<-1.?-1.:x,irts)
        
    if length(irts)==0
        Float64[]
    else
        fromcanonical(f,irts)
    end
end


##TODO: allow using routines for complex domains below
function Base.maximum(f::IFun)
    pts=[f.domain.a,f.domain.b,roots(diff(f))]
    maximum(f[pts])
end

function Base.minimum(f::IFun)
    pts=[f.domain.a,f.domain.b,roots(diff(f))]
    minimum(f[pts])
end

function Base.indmax(f::IFun)
    pts=[f.domain.a,f.domain.b,roots(diff(f))]
    pts[indmax(f[pts])]
end

function Base.indmin(f::IFun)
    pts=[f.domain.a,f.domain.b,roots(diff(f))]
    pts[indmin(f[pts])]
end



