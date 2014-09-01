export TensorFun



type TensorFun{F<:IFun,D<:IntervalDomain}
    coefficients::Vector{F}     # coefficients are in x
    domainy::D
end

function TensorFun{T<:Number}(cfs::Array{T},dx,dy)
    ret=Array(IFun{T,typeof(dx)},size(cfs,2))
    for k=1:size(cfs,2)
        ret[k]=IFun(chop!(cfs[:,k],10eps()),dx)
    end
    TensorFun(ret,dy)
end

TensorFun(f::Fun2D)=TensorFun(coefficients(f),domain(f,1),domain(f,2))

TensorFun(f::Function,d1...)=TensorFun(Fun2D(f,d1...))

Base.size(f::TensorFun,k::Integer)=k==1?mapreduce(length,max,f.coefficients):length(f.coefficients)
Base.size(f::TensorFun)=(size(f,1),size(f,2))


function funlist2coefficients{T<:Number,D<:IntervalDomain}(f::Vector{IFun{T,D}})
    A=zeros(T,mapreduce(length,max,f),length(f))
    for k=1:length(f)
        A[1:length(f[k]),k]=f[k].coefficients
    end
end

coefficients(f::TensorFun)=funlist2coefficients(f.coefficients)

function coefficients(f::TensorFun,ox::Integer,oy::Integer)
    m=size(f,1)
    A=[pad!(coefficients(fx,ox),m) for fx in f.coefficients]
    B=hcat(A...)::Array{Float64,2}
    for k=1:size(B,1)
        B[k,:]=ultraconversion(vec(B[k,:]),oy)
    end
    
    B
end


domain(f::TensorFun,k::Integer)=k==1?domain(f.coefficients[1]):f.domainy




Fun2D(f::TensorFun)=Fun2D(f.coefficients.',domain(f,2))

.'(f::TensorFun)=TensorFun(coefficients(f).',domain(f,2),domain(f,1))
