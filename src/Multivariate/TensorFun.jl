export TensorFun



type TensorFun{F<:IFun,D<:IntervalDomain}
    coefficients::Vector{F}##coefficients are in x
    domainy::D
end

function TensorFun(f::Fun2D)
    cfs=coefficients(f)
    dx=domain(f,1)
    dy=domain(f,2)
    ret=Array(typeof(f.A[1]),size(cfs,2))
    for k=1:size(cfs,2)
        ret[k]=IFun(cfs[:,k],dx)
    end
    TensorFun(ret,dy)
end


TensorFun(f::Function,d1...)=TensorFun(Fun2D(f,d1...))

Base.size(f::TensorFun,k::Integer)=k==1?mapreduce(length,max,f.coefficients):length(f.coefficients)
Base.size(f::TensorFun)=(size(f,1),size(f,2))

function coefficients{T<:Number}(f::TensorFun{IFun{T}})
    A=zeros(T,size(f)...)
    for k=1:length(f.coefficients)
        A[1:length(f.coefficients[k]),k]=f.coefficients[k].coefficients
    end
end 

domain(f::TensorFun,k::Integer)=k==1?domain(f.coefficients[1]):f.domainy




Fun2D(f::TensorFun)=Fun2D(f.coefficients,domain(f,2))