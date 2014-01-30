include("ShiftVector.jl")

##TODO: unify tolerances



type FFun{T<:Number,D<:PeriodicDomain} <: AbstractFun
    coefficients::ShiftVector{T}
    domain::D
end



FFun(f::Function,n::Integer)=FFun(f,PeriodicInterval(),n)
FFun(f::Function,d::Domain,n::Integer)=FFun(svfft(f(points(d,n))),d)
FFun(f::Function,d::Vector,n::Integer)=FFun(f,apply(PeriodicInterval,d),n)
FFun(cfs::ShiftVector)=FFun(cfs,PeriodicInterval())
FFun(cfs::ShiftVector,d::Vector)=FFun(cfs,apply(PeriodicInterval,d))
FFun(f::Function)=FFun(f,PeriodicInterval())
FFun(f::Function,d::Vector)=FFun(f,apply(PeriodicInterval,d))

function FFun(f::Function, d::Domain)
    #TODO: reuse function values

    tol = 200*eps();

    oldcf = FFun(f,d,2);

    for logn = 2:20
        cf = FFun(f, d, 2^logn);
        
        if max(
            abs(last(cf.coefficients)),
            abs(cf.coefficients.vector[end-1]),
            abs(cf.coefficients.vector[1]),
            abs(cf.coefficients.vector[2])
##TODO: compare coefficients            
#             ,
#             max(abs(cf.coefficients[1:2^(logn-1) + 1] - oldcf.coefficients))
        ) < tol
            return cf;
        end
        
        oldcf = cf;
    end
    
    warn("Maximum length reached");
    
    cf
end





##Evaluation


evaluate(f::FFun,x)=horner(f.coefficients,tocanonical(f,x))


##TODO: fast list routine
function horner(v::ShiftVector,x)
    ret = 0.;
    ei = exp(1.im.*x);
    ein = exp(-1.im.*x);    
    
    p = 1.;
    for k = 0:lastindex(v)
        ret += v[k]*p;
        p .*= ei;
    end
    p=ein;
    for k = -1:-1:firstindex(v)
        ret+= v[k]*p;
        p .*= ein;
    end
    
    ret
end

Base.getindex(f::FFun,x)=evaluate(f,x)

##Data routines  TODO: Unify with IFun
values(f::FFun)=isvfft(f.coefficients)
points(f::FFun)=points(f.domain,length(f))
Base.length(f::FFun)=length(f.coefficients)


## Manipulate length

pad(f::FFun,r::Range1)=FFun(pad(f.coefficients,r),f.domain)


## Addition

## Addition and multiplication




for op = (:+,:-)
    @eval begin
        function ($op)(f::FFun,g::FFun)
            @assert f.domain == g.domain

            fi = min(firstindex(f.coefficients),firstindex(g.coefficients))        
            li = max(lastindex(f.coefficients),lastindex(g.coefficients))
            f2 = pad(f,fi:li);
            g2 = pad(g,fi:li);
            
            FFun(($op)(f2.coefficients,g2.coefficients),f.domain)
        end

        function ($op)(f::FFun,c::Number)
            f2 = deepcopy(f);
            
            f2.coefficients[0] = ($op)(f2.coefficients[0],c);
            
            f2
        end
    end
end 



function .*(f::FFun,g::FFun)
    @assert f.domain == g.domain
    #TODO Coefficient space version
    fi = firstindex(f.coefficients)+firstindex(g.coefficients)
    li = lastindex(f.coefficients)+lastindex(g.coefficients)
    
    ##Todo: Fix following hack to ensure FFT index range
    li = max(-fi,li)
    fi = min(fi,-li)    
    
    
    f2 = pad(f,fi:li)
    g2 = pad(g,fi:li)
    

    FFun(svfft(values(f2).*values(g2)),f.domain)
end




for op = (:*,:.*,:./,:/)
    @eval ($op)(f::FFun,c::Number) = FFun(($op)(f.coefficients,c),f.domain)
end 

-(f::FFun)=FFun(-f.coefficients,f.domain)
-(c::Number,f::FFun)=-(f-c)


for op = (:*,:.*,:+)
    @eval ($op)(c::Number,f::FFun)=($op)(f,c)
end


##TODO: Division by FFun

## Norm

import Base.norm


##TODO: this is mapped L^2[-π,π] norm, 
norm(f::FFun)=norm(f.coefficients.vector)


##TODO: Implement coefficient space real and imag

import Base.imag, Base.real

for op = (:real,:imag) 
    @eval ($op)(f::FFun) = IFun(svfft(($op)(values(f))),f.domain)
end



## Integration

Base.first(cf::FFun)=(-1.)^(cf.coefficients.index+1).*foldr(-,cf.coefficients.vector)

function Base.cumsum(f::FFun)
    cf = integrate(f)
    cf - first(cf)
end




##TODO: Root finding

##TODO: Sampling

