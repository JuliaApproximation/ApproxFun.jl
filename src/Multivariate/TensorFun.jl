export TensorFun



immutable TensorFun{T<:Union(Float64,Complex{Float64}),S<:FunctionSpace,V<:FunctionSpace}<:MultivariateFun
    coefficients::Vector{Fun{T,S}}     # coefficients are in x
    spacey::V
end

for T in (:Float64,:(Complex{Float64}))
    @eval TensorFun{F<:Fun{$T}}(M::Vector{F},dy::FunctionSpace)=TensorFun{$T,typeof(M[1].space),typeof(dy)}(Fun{$T,typeof(M[1].space)}[Mk for Mk in M],dy)
end

function TensorFun{T<:Number,S<:FunctionSpace}(cfs::Matrix{T},dx::S,dy::FunctionSpace)
    ret=Array(Fun{T,S},size(cfs,2))
    for k=1:size(cfs,2)
        ret[k]=chop!(Fun(cfs[:,k],dx),10eps())
    end
    TensorFun{T,S,typeof(dy)}(ret,dy)
end

TensorFun(cfs::Array,d::TensorSpace)=TensorFun(cfs,d[1],d[2])
TensorFun(cfs::Array,d::ProductDomain)=TensorFun(cfs,d[1],d[2])
TensorFun(f::Function,dy::Domain)=error("This function is only implemented to avoid ambiguity, do not call.")
TensorFun(f,dy::Domain)=TensorFun(f,Space(dy))
TensorFun(f,dx::Domain,dy::Domain)=TensorFun(f,Space(dx),Space(dy))
TensorFun(f::Fun2D)=TensorFun(coefficients(f),space(f,1),space(f,2))

TensorFun(f::Function,d1...)=TensorFun(Fun2D(f,d1...))

Base.size(f::TensorFun,k::Integer)=k==1?mapreduce(length,max,f.coefficients):length(f.coefficients)
Base.size(f::TensorFun)=(size(f,1),size(f,2))

for T in (:Float64,:(Complex{Float64}))
    @eval begin
        function funlist2coefficients{F<:Fun{$T}}(f::Vector{F})
            A=zeros($T,mapreduce(length,max,f),length(f))
            for k=1:length(f)
                A[1:length(f[k]),k]=f[k].coefficients
            end
            A
        end
    end
end


coefficients(f::TensorFun)=funlist2coefficients(f.coefficients)

function coefficients{T<:Number}(f::TensorFun{T},ox::FunctionSpace,oy::FunctionSpace)
    m=size(f,1)
    A=[pad!(coefficients(fx,ox),m) for fx in f.coefficients]
    B=hcat(A...)::Array{T,2}
    for k=1:size(B,1)
        B[k,:]=spaceconversion(vec(B[k,:]),space(f,2),oy)
    end
    
    B
end

function values(f::TensorFun)
    n=size(f,1)
    M=hcat(map(fk->values(pad(fk,n)),f.coefficients)...)
    for k=1:size(f,1)
        M[k,:]=values(Fun(vec(M[k,:]),space(f,2)))
    end
    M     
end




points(f::TensorFun,k)=points(space(f,k),size(f,k))

space(f::TensorFun,k::Integer)=k==1?space(f.coefficients[1]):f.spacey
domain(f::TensorFun,k::Integer)=domain(space(f,k))



evaluate{T}(f::TensorFun{T},x::Real,::Colon)=Fun(T[fc[x] for fc in f.coefficients],space(f,2))
evaluate(f::TensorFun,x::Real,y::Real)=evaluate(f,x,:)[y]
evaluate(f::TensorFun,x::Colon,y::Real)=evaluate(f.',y,:)
evaluate(f::TensorFun,xx::Vector,yy::Vector)=hcat([evaluate(f,x,:)[[yy]] for x in xx]...).'
evaluate(f::TensorFun,x::Range,y::Range)=evaluate(f,[x],[y])


*(c::Number,f::TensorFun)=TensorFun(c*f.coefficients,f.spacey)
*(f::TensorFun,c::Number)=c*f



chop(f::TensorFun,es...)=TensorFun(map(g->chop(g,es...),f.coefficients),f.spacey)


##TODO: following assumes f is never changed....maybe should be deepcopy?
function +(f::TensorFun,c::Number)
    cfs=copy(f.coefficients)
    cfs[1]+=c
    TensorFun(cfs,f.spacey)
end
+(c::Number,f::TensorFun)=f+c
-(f::TensorFun,c::Number)=f+(-c)
-(c::Number,f::TensorFun)=c+(-f)


function +(f::TensorFun,g::TensorFun)
    if size(f,2) >= size(g,2)
        @assert f.spacey==g.spacey
        cfs = copy(f.coefficients)
        for k=1:size(g,2)
            cfs[k]+=g.coefficients[k]
        end
        
        TensorFun(cfs,f.spacey)
    else
        g+f
    end
end

-(f::TensorFun)=(-1)*f
-(f::TensorFun,g::TensorFun)=f+(-g)


Fun2D(f::TensorFun)=Fun2D(f.coefficients,space(f,2))

function Base.diff(f::TensorFun,j::Integer)
    if j==1
        TensorFun(map(diff,f.coefficients),f.spacey)
    else
        diff(f.',1).'
    end
end


Base.transpose(f::TensorFun)=TensorFun(coefficients(f).',space(f,2),space(f,1))





#TODO: adaptive
#TODO: assumes Chebyshev
for op in (:(Base.sin),:(Base.cos))
    @eval ($op)(f::TensorFun)=TensorFun(chebyshevtransform($op(values(f))),domain(f))
end

#TODO: assumes real basis
for op = (:(Base.real),:(Base.imag),:(Base.conj)) 
    @eval ($op)(f::TensorFun) = TensorFun(map($op,f.coefficients),f.spacey)
end

