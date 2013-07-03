## Helper routines



alternating_vector(n::Integer) = 2*mod([1:n],2)-1;

function chebyshev_transform(x::Vector)
    ret = FFTW.r2r(x, FFTW.REDFT00);
    ret[1] *= .5;
    ret[end] *= .5;    
    ret.*alternating_vector(length(ret))/(length(ret)-1)
end

function ichebyshev_transform(x::Vector)
    x[1] *= 2.;
    x[end] *= 2.;
    
    ret = chebyshev_transform(x);
    
    x[1] *= .5;
    x[end] *= .5;
    
    ret[1] *= 2.;
    ret[end] *= 2.;
    
    flipud(ret.*alternating_vector(length(ret)).*(length(x) - 1).*.5)
end

points(n::Integer)= cos(π*[n-1:-1:0]/(n-1))

points(d::Domain,n::Integer) = from_uinterval(d,points(n))



##  Constructors

abstract Fun



type IFun{T<:Number} <: Fun
    coefficients::Vector{T}
    domain::Domain
end



function IFun(f::Function,n::Integer)
    IFun(f,Interval(-1,1),n)
end

function IFun(f::Function,d::Domain,n::Integer)
    IFun(chebyshev_transform(f(points(d,n))),d)
end

function IFun(f::Function,d::Vector,n::Integer)
    IFun(f,apply(Interval,d),n)
end

function IFun(cfs::Vector)
	IFun(cfs,Interval(-1,1))
end

function IFun(cfs::Vector,d::Vector)
	IFun(cfs,apply(Interval,d))
end


function IFun(f::Function)
    IFun(f,Interval(-1,1))
end

function IFun(f::Function,d::Vector)
    IFun(f,apply(Interval,d))
end

function IFun(f::Function, d::Domain)
    #reuse function values

    tol = 200*eps();

    for logn = 2:20
        cf = IFun(f, d, 2^logn);
        
        if max(abs(cf.coefficients[end]),abs(cf.coefficients[end-1])) < tol
            chop!(cf,tol);
            return cf;
        end
    end
    
    warn("Maximum length reached");
    
    cf
end



function evaluate(f::IFun,x)
    evaluate(f.coefficients,to_uinterval(f.domain,x))
end


##TODO: Speed up vector method

function evaluate(v::Vector,xl::Vector)
    map(x->evaluate(v,x),xl)
end

function evaluate(v::Vector,x::Real)
    @assert (-1 <= x) && (x <= 1)

    unp = 0.;
    un = v[end];
    n = length(v);
    for k = n-1:-1:2
        uk = 2.*x.*un - unp + v[k];
        unp = un;
        un = uk;
    end

    uk = 2.*x.*un - unp + 2*v[1];
    .5*(uk -unp)
end



function values(f::IFun)
   ichebyshev_transform(f.coefficients) 
end

function points(f::IFun)
    points(f.domain,length(f))
end


function Base.length(f::IFun)
    length(f.coefficients)
end


## Matipulate length


function pad!(f::IFun,n::Integer)
	if (n > length(f))
		append!(f.coefficients,zeros(n - length(f)));
	else
		resize!(f.coefficients,n);
	end
end


function pad(f::IFun,n::Integer)
	if (n > length(f))
		IFun([f.coefficients,zeros(n - length(f))],f.domain);
	else
		IFun(f.coefficients[1:n],f.domain);
	end
end


function chop!(f::IFun,tol::Real)
    @assert tol > 0

    for k=[length(f):-1:1]
        if abs(f.coefficients[k]) > tol
            resize!(f.coefficients,k);
            return;
        end
    end
end


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
    
    IFun(chebyshev_transform(values(f2).*values(g2)),f.domain)
end




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
            
            f2 = IFun(chebyshev_transform(c/values(f2)),f.domain);
            
            if max(abs(f2.coefficients[end-1:end])) > 10*eps()
                warn("Division has not converged, may be inaccurate")
            end
            
            f2
        end
    end
end

./(f::IFun,g::IFun)=f.*(1./g)



## Norm


## Mapped functions

import Base.imag, Base.real

for op = (:real,:imag) 
    @eval ($op)(f::IFun) = IFun(($op)(f.coefficients),f.domain)
end


## Differentiation and integration


## Start of support for UFun

# diff from T -> U
ultradiff(v::Vector)=[1:length(v)-1].*v[2:end]

#int from U ->T
ultraint(v::Vector)=[0,v./[1:length(v)]]



# Convert from U -> T
function ultraiconv(v::Vector)
    n = length(v);
    w = zeros(n);
    
    w[end] = 2v[end];
    w[end-1] = 2v[end-1];
    
    for k = n-2:-1:2
        w[k] = 2*(v[k] + .5w[k+2]);
    end
    
    w[1] = v[1] + .5w[3];
    
    w
end


# Convert U to T
function ultraconv(v::Vector)
    n = length(v);
    w = zeros(n);
    
    if n == 1
        w[1] = v[1];
    elseif n == 2
        w[1] = v[1];
        w[2] = .5v[2];
    else
        w[1] = v[1] - .5v[3];        
    
        for k = 2:n-2
            w[k] = .5*(v[k] - v[k+2]);
        end
    
        w[n-1] = .5v[n-1];
        w[n] = .5v[n];        
    end
    
    w
end


# diff T -> U, then convert U -> T
function Base.diff(f::IFun)

    # Will need to change code for other domains
    @assert typeof(f.domain) <: Interval
    
    to_uintervalD(f.domain,0)*IFun(ultraiconv(ultradiff(f.coefficients)),f.domain)
end

function Base.cumsum(f::IFun)
    # Will need to change code for other domains
    @assert typeof(f.domain) <: Interval
    
    from_uintervalD(f.domain,0)*IFun(ultraint(ultraconv(f.coefficients)),f.domain)    
end


==(f::IFun,g::IFun) =  (f.coefficients == g.coefficients && f.domain == g.domain)



## Root finding


## Plotting


# function Winston.plot(f::IFun{Complex})
#     plot(points(f),values(real(f)),points(f),values(imag(f)))
# end


Winston.plot(f::IFun) = plot(points(f),values(real(f)),points(f),values(imag(f)))







