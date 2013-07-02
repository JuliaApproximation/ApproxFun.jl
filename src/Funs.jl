module Funs

export IFun,Interval,evaluate,values,points,chebyshev_transform

alternating_vector(n::Integer) = 2*mod([1:n],2)-1;

function chebyshev_transform(x::Vector)
    ret = FFTW.r2r(x, FFTW.REDFT00);
    ret[1] *= .5;
    ret[end] *= .5;    
    ret.*alternating_vector(length(ret))/(length(ret)-1)
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


function evaluate(f,x)
    f.coefficients
end



function values(f)
    f.domain
end

end #module
