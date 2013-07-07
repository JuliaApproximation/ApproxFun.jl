

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


function evaluate(f::FFun,x)
    0
end
Base.getindex(f::FFun,x)=evaluate(f,x)


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
