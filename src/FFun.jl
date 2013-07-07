

periodic_interval = PeriodicInterval(-1.π,1.π)

type FFun{T<:Number} <: AbstractFun
    coefficients::ShiftVector{T}
    domain::Domain
end


function FFun(f::Function,n::Integer)
    FFun(f,periodic_interval,n)
end

function FFun(f::Function,d::Domain,n::Integer)
    FFun(svfft(f(points(d,n))),d)
end

function FFun(f::Function,d::Vector,n::Integer)
    FFun(f,apply(PeriodicInterval,d),n)
end

function FFun(cfs::ShiftVector)
	FFun(cfs,periodic_interval)
end

function FFun(cfs::ShiftVector,d::Vector)
	FFun(cfs,apply(PeriodicInterval,d))
end


function FFun(f::Function)
    FFun(f,periodic_interval)
end

function FFun(f::Function,d::Vector)
    FFun(f,apply(PeriodicInterval,d))
end


##Evaluation


evaluate(f::FFun,x)=horner(f.coefficients,tocanonical(f,x))

function horner(v::ShiftVector,x)
    ret = 0.;
    ei = exp(1.im.*x);
    ein = exp(-1.im.*x);    
    
    p = 1.;
    for k = 0:lastindex(v)
        ret += v[k]*p;
        p *= ei;
    end
    p=ein;
    for k = -1:-1:firstindex(v)
        ret+= v[k]*p;
        p *= ein;
    end
    
    ret
end

Base.getindex(f::FFun,x)=evaluate(f,x)

##Data routines  TODO: Unify with IFun
values(f::FFun)=0
points(f::IFun)=points(f.domain,length(f))
Base.length(f::IFun)=length(f.coefficients)


##Plotting

function plot(f::FFun{Float64}) 
    p = FramedPlot();
    pts = points(pad(f,2length(f)));
    vals =values(pad(f,2length(f)));
    add(p,Curve(pts,vals,"color","blue"));   
    Winston.display(p);
    p
end

function plot(f::FFun{Complex{Float64}}) 
    p = FramedPlot();
    pts = points(pad(f,2length(f)));
    vals =values(pad(f,2length(f)));
    add(p,Curve(pts,real(vals),"color","blue"));
    add(p,Curve(pts,imag(vals),"color","red"));    
    Winston.display(p);
    p
end
