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