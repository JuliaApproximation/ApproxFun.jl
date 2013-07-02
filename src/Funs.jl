module Funs
    using Base

export IFun,Interval,evaluate,values,points,chebyshev_transform

alternating_vector(n::Integer) = 2*mod([1:n],2)-1;

function chebyshev_transform(x::Vector)
    ret = FFTW.r2r(x, FFTW.REDFT00);
    ret[1] *= .5;
    ret[end] *= .5;    
    ret.*alternating_vector(length(ret))/(length(ret)-1)
end

function ichebyshev_transform(x::Vector)
    x[1] *= 2;
    x[end] *= 2;
    
    ret = chebyshev_transform(x);
    
    x[1] *= .5;
    x[end] *= .5;
    
    ret[1] *= 2;
    ret[end] *= 2;
    
    flipud(ret.*alternating_vector(length(ret)).*(length(x) - 1).*.5)
end

points(n)= cos(π*[n-1:-1:0]/(n-1))

# points(d::Interval,n) = (d == [-1,1]) ? points(n) : points(n)


type IFun
    coefficients::Vector
    domain
end



function IFun(f::Function,n::Integer)
    IFun(chebyshev_transform(f(points(n))),[-1,1])
end


function evaluate(f::IFun,x)
    unp = 0;
    un = f.coefficients[end];
    n = length(f);
    for k = n-1:-1:2
        uk = 2.*x.*un - unp + f.coefficients[k];
        unp = un;
        un = uk;
    end

    uk = 2.*x.*un - unp + 2*f.coefficients[1];
    .5*(uk -unp)
end



function values(f::IFun)
   ichebyshev_transform(f.coefficients) 
end

function Base.length(f::IFun)
    length(f.coefficients)
end


end #module
