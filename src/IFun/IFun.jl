include("bary.jl")
include("clenshaw.jl")
include("ultraspherical.jl")


export plan_chebyshevtransform

## transforms


function negateeven!(x)
    for k =2:2:length(x)
        x[k] *= -1.
    end
    
    x
end

plan_chebyshevtransform(x)=FFTW.plan_r2r(x, FFTW.REDFT00)
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




IFun(f::Function,n::Integer)=IFun(f,Interval(),n)
IFun(f::Function,d::Domain,n::Integer)=IFun(chebyshevtransform(f(points(d,n))),d)
IFun(f::Function,d::Vector,n::Integer)=IFun(f,apply(Interval,d),n)
IFun(cfs::Vector)=IFun(1.0*cfs,Interval())
IFun(cfs::Vector,d::Vector)=IFun(1.0*cfs,apply(Interval,d))
IFun(cfs::Vector,d::IntervalDomain)=IFun(1.0*cfs,d)
IFun(f::Function)=IFun(f,Interval())
IFun(f::Function,d::Vector)=IFun(f,apply(Interval,d))

IFun(f::IFun,d)=IFun(f.coefficients,d)

## Adaptive constructors

function randomIFun(f::Function,d::Domain)
    @assert d == Interval()

    #TODO: implement other domains
    
    IFun(chebyshevtransform(randomadaptivebary(f)),d)
end

function zerocfsIFun(f::Function,d::Domain)
    #reuse function values

    tol = 200*eps();

    oldcf = IFun(f,d,2 + 1);

    for logn = 2:20
        cf = IFun(f, d, 2^logn + 1);
        
        if max(abs(cf.coefficients[end]),abs(cf.coefficients[end-1]),maximum(abs(cf.coefficients[1:2^(logn-1) + 1] - oldcf.coefficients))) < tol
            chop!(cf,10eps());
            return cf;
        end
        
        oldcf = cf;
    end
    
    warn("Maximum length reached");
    
    oldcf
end


function IFun(f::Function, d::Domain; method="zerocoefficients")
    if method == "zerocoefficients"
        zerocfsIFun(f,d)
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



## Matipulate length


pad!(f::IFun,n::Integer)=pad!(f.coefficients,n)
pad(f::IFun,n::Integer)=IFun(pad(f.coefficients,n),f.domain)





function chop!(c::Vector,tol::Real)
    @assert tol > 0

    for k=[length(c):-1:1]
        if abs(c[k]) > tol
            resize!(c,k);
            return c;
        end
    end
    
    []
end

chop!(f::IFun,tol::Real)=(chop!(f.coefficients,tol); f)
Base.chop(f,tol)=chop!(deepcopy(f),tol)


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
    
    IFun(chebyshevtransform(values(f2).*values(g2)),f.domain)
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




# division by fun 

for op = (:./,:/)
    @eval begin
        function ($op)(c::Number,f::IFun)
            #TODO choose the length
        
            f2 = pad(f,2*length(f));
            
            f2 = IFun(chebyshevtransform(c/values(f2)),f.domain);
            
            if maximum(abs(f2.coefficients[end-1:end])) > 10*eps()
                warn("Division has not converged, may be inaccurate")
            end
            
            f2
        end
    end
end

./(f::IFun,g::IFun)=f.*(1./g)


##Coefficient space operators


function multiplybyx(f::IFun)
    a = f.domain.a
    b = f.domain.b
    g = IFun([0,1,.5*ones(length(f)-1)].*[0,f.coefficients]+[.5*f.coefficients[2:end],0,0],f.domain) #Gives multiplybyx on unit interval
    (b-a)/2*g + (b+a)/2
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
    c=chop(cin,10eps());
    if c == [] || length(c) == 1
        return []
    elseif length(c) == 2
        return [-c[1]/c[2]]
    else 
        n=length(c)-1;
        
        I = [ones(Int64,n),2:n-1,2:n];
        J=[1:n,3:n,1:n-1];
        V = [-c[end-1]/(2c[end]),.5-c[end-2]/(2c[end]),-c[end-3:-1:1]/(2c[end]),.5*ones(n-2),.5*ones(n-2),1];
        C=sparse(I,J,V);
        A=zeros(n,n);
        A[1:end,1:end]=C[1:end,1:end];
        
        Λ,V=eig(A);
        return Λ    
    end
end


function complexroots(f::IFun)
    fromcanonical(f,complexroots(f.coefficients))
end

function roots(f::IFun)
    complexroots(f)
end

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


    



## Plotting


# function Winston.plot(f::IFun{Complex})
#     plot(points(f),values(real(f)),points(f),values(imag(f)))
# end


using Winston

function plot(f::IFun{Float64}) 
    p = Winston.FramedPlot();
    pts = points(pad(f,2length(f)));
    vals =values(pad(f,2length(f)));
    add(p,Curve(pts,vals,"color","blue"));   
    Winston.display(p);
#    p
end

function plot(f::IFun{Complex{Float64}}) 
    p = Winston.FramedPlot();
    pts = points(pad(f,2length(f)));
    vals =values(pad(f,2length(f)));
    add(p,Curve(pts,real(vals),"color","blue"));
    add(p,Curve(pts,imag(vals),"color","red"));    
    Winston.display(p);
#    p
end




## Array routines

function values{T,D}(p::Array{IFun{T,D},1})
    n = maximum(map(length,p))
    ret = Array(T,length(p),n)
    for i = 1:length(p)
        ret[i,:] = values(pad(p[i],n));
    end
    ret
end

function values{T,D}(p::Array{IFun{T,D},2})
    @assert size(p)[1] == 1

    n = maximum(map(length,p))
    ret = Array(T,n,length(p))
    for i = 1:length(p)
        ret[:,i] = values(pad(p[i],n));
    end
    ret
end

function coefficients{T,D}(p::Array{IFun{T,D},1})
    n = maximum(map(length,p))
    ret = Array(T,length(p),n)
    for i = 1:length(p)
        ret[i,:] = pad(p[i],n).coefficients
    end
    ret
end




function *{T,D}(A::Array{Float64,2},p::Array{IFun{T,D},1})
    cfs=(A*coefficients(p))'
    ret = Array(IFun{T,D},size(A)[1])
    for i = 1:size(A)[1]
        ret[i] = IFun(cfs[:,i],p[i].domain)
    end
    ret
end

