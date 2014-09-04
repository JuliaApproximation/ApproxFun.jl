export TensorFun



type TensorFun{F<:IFun,D<:IntervalDomain}<:MultivariateFun
    coefficients::Vector{F}     # coefficients are in x
    domainy::D
end

function TensorFun{T<:Number}(cfs::Array{T},dx,dy)
    ret=Array(IFun{T,typeof(dx)},size(cfs,2))
    for k=1:size(cfs,2)
        ret[k]=chop!(IFun(cfs[:,k],dx),10eps())
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
    A
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

function values(f::TensorFun)
    m,n=size(f)
    p=plan_chebyshevtransform(pad(f.coefficients[1].coefficients,m))
    p2=plan_chebyshevtransform(pad(f.coefficients[1].coefficients,n))
    A=Array(Float64,m,n)
    for k=1:n
        A[:,k]=ichebyshevtransform(pad(f.coefficients[k].coefficients,m),p)
    end
    
    for k=1:m
        A[k,:]=ichebyshevtransform(vec(A[k,:]),p2)
    end
    
    
    A
end

points(f::TensorFun,k)=points(domain(f,k),size(f,k))

domain(f::TensorFun,k::Integer)=k==1?domain(f.coefficients[1]):f.domainy
domain(LL::TensorFun)=domain(LL,1)⊗domain(LL,2)

evaluate(f::TensorFun,x::Real,::Colon)=IFun([fc[x] for fc in f.coefficients],f.domainy)
evaluate(f::TensorFun,x::Real,y::Real)=evaluate(f,x,:)[y]
evaluate(f::TensorFun,x::Colon,y::Real)=evaluate(f.',y,:)
evaluate(f::TensorFun,xx::Vector,yy::Vector)=hcat([evaluate(f,x,:)[[yy]] for x in xx]...).'
evaluate(f::TensorFun,x::Range,y::Range)=evaluate(f,[x],[y])


*(c::Number,f::TensorFun)=TensorFun(c*f.coefficients,f.domainy)
*(f::TensorFun,c::Number)=c*f



chop(f::TensorFun,es...)=TensorFun(map(g->chop(g,es...),f.coefficients),f.domainy)


##TODO: following assumes f is never changed....maybe should be deepcopy?
function +(f::TensorFun,c::Number)
    cfs=copy(f.coefficients)
    cfs[1]+=c
    TensorFun(cfs,f.domainy)
end
+(c::Number,f::TensorFun)=f+c
-(f::TensorFun,c::Number)=f+(-c)
-(c::Number,f::TensorFun)=c+(-f)


function +(f::TensorFun,g::TensorFun)
    if size(f,2) >= size(g,2)
        @assert f.domainy==g.domainy
        cfs = copy(f.coefficients)
        for k=1:size(g,2)
            cfs[k]+=g.coefficients[k]
        end
        
        TensorFun(cfs,f.domainy)
    else
        g+f
    end
end

-(f::TensorFun)=(-1)*f
-(f::TensorFun,g::TensorFun)=f+(-g)


Fun2D(f::TensorFun)=Fun2D(f.coefficients.',domain(f,2))

.'(f::TensorFun)=TensorFun(coefficients(f).',domain(f,2),domain(f,1))



