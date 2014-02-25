module ApproxFun
    using Base

export Fun,IFun,FFun,Interval,evaluate,values,points,chebyshevtransform
export pad!,pad,sample,chop!,complexroots,roots,plot,svfft
export multiplybyx,IntervalDomain,fasttimes

##Testing
export bisectioninv, clenshaw
export coefficients, integrate



abstract AbstractFun


## Helper routines
alternatingvector(n::Integer) = 2*mod([1:n],2)-1

function pad!(f::Vector,n::Integer)
	if (n > length(f))
		append!(f,zeros(n - length(f)))
	else
		resize!(f,n)
	end
end


function pad(f::Vector,n::Integer)
	if (n > length(f))
        [f,zeros(n - length(f))]
	else
        f[1:n]
	end
end

function pad{T}(A::Array{T},n::Integer,m::Integer)
	if n <= size(A,1) && m <= size(A,2)
        A[1:n,1:m]
	else
        ret = zeros(T,n,m)
        
        if n <= size(A,1)
            ret[1:n,1:size(A,2)] = A[1:n,:]
        elseif m <= size(A,2)
            ret[1:size(A,1),1:m] = A[:,1:m]
        else
            ret[1:size(A,1),1:size(A,2)]=A
        end
        
        ret
	end
end




#TODO:padleft!

function padleft(f::Vector,n::Integer)
	if (n > length(f))
        [zeros(n - length(f)),f]
	else
        f[end-n+1:end]
	end
end



include("Domains/Domain.jl")



##Chebyshev Routines
include("IFun/IFun.jl")


##Fourier Routines
include("FFun/FFun.jl")

# Canonical domains
include("Domains/Interval.jl")
include("Domains/PeriodicInterval.jl")


include("Operators/Operator.jl")


include("Multivariate/VectorFun.jl")
include("Multivariate/Fun2D.jl")

include("sample.jl")


##Following routine decides
# whether input is IFun or FFun
Fun(x)=IFun(x)
Fun(x,d::IntervalDomain)=IFun(x,d)
Fun(x,d::PeriodicDomain)=FFun(x,d)
Fun(x,d)=IFun(x,d)
Fun(x,d::IntervalDomain,n::Integer)=IFun(x,d,n)
Fun(x,d::PeriodicDomain,n::Integer)=FFun(x,d,n)
Fun(x,d::Vector,n::Integer)=IFun(x,d,n)



## Other domains

include("Domains/Ray.jl")
include("Domains/Line.jl")
include("Domains/Circle.jl")

## Further extra features

include("Singularities/SingFun.jl")

include("PDE/pdesolve.jl")

include("Plotting/IFunPlot.jl")





end #module


