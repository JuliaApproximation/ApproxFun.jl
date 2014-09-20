export TensorFun,ProductFun


abstract AbstractProductFun{S,V,T} <: MultivariateFun

for TY in (:TensorFun,:ProductFun)
    @eval immutable $TY{S<:FunctionSpace,V<:FunctionSpace,T<:Union(Float64,Complex{Float64})}<:AbstractProductFun{S,V,T}
        coefficients::Vector{Fun{S,T}}     # coefficients are in x
        spacey::V
    end
end
for TY in (:TensorFun,:ProductFun), T in (:Float64,:(Complex{Float64}))
    @eval $TY{S}(M::Vector{Fun{S,$T}},dy::FunctionSpace)=$TY{$T,typeof(M[1].space),typeof(dy)}(Fun{typeof(M[1].space),$T}[Mk for Mk in M],dy)
end


function TensorFun{T<:Number,S<:FunctionSpace}(cfs::Matrix{T},dx::S,dy::FunctionSpace)
    ret=Array(Fun{S,T},size(cfs,2))
    for k=1:size(cfs,2)
        ret[k]=chop!(Fun(cfs[:,k],dx),10eps())
    end
    TensorFun{S,typeof(dy),T}(ret,dy)
end

function ProductFun{T<:Number,S<:FunctionSpace}(cfs::Matrix{T},dx::Vector{S},dy::FunctionSpace)
    @assert size(cfs,2)==length(dx)
    
    ret=Array(Fun{S,T},size(cfs,2))
    for k=1:size(cfs,2)
        ret[k]=chop!(Fun(cfs[:,k],dx[k]),10eps())
    end
    TensorFun{S,typeof(dy),T}(ret,dy)
end

TensorFun(cfs::Array,d::TensorSpace)=TensorFun(cfs,d[1],d[2])
TensorFun(cfs::Array,d::ProductDomain)=TensorFun(cfs,d[1],d[2])
TensorFun(f::Function,dy::Domain)=error("This function is only implemented to avoid ambiguity, do not call.")
TensorFun(f,dy::Domain)=TensorFun(f,Space(dy))
TensorFun(f,dx::Domain,dy::Domain)=TensorFun(f,Space(dx),Space(dy))
TensorFun(f::LowRankFun)=TensorFun(coefficients(f),space(f,1),space(f,2))

TensorFun(f::Function,d1...)=TensorFun(LowRankFun(f,d1...))
function ProductFun{ST<:FunctionSpace}(f::Function,S::Vector{ST},T::FunctionSpace,N::Integer)
    M=length(S)     # We assume the list of spaces is the same length
    xx=points(S[1],N)
    tt=points(T,M)
    V=Float64[f(x,t) for x=xx, t=tt]
    ProductFun(transform(S,T,V),S,T)
end

# For specifying spaces by anonymous function
ProductFun(f::Function,SF::Function,T::FunctionSpace,N::Integer,M::Integer)=ProductFun(f,typeof(SF(1))[SF(k) for k=1:M],T,N)

Base.size(f::AbstractProductFun,k::Integer)=k==1?mapreduce(length,max,f.coefficients):length(f.coefficients)
Base.size(f::AbstractProductFun)=(size(f,1),size(f,2))

for T in (:Float64,:(Complex{Float64}))
    @eval begin
        function funlist2coefficients{S}(f::Vector{Fun{S,$T}})
            A=zeros($T,mapreduce(length,max,f),length(f))
            for k=1:length(f)
                A[1:length(f[k]),k]=f[k].coefficients
            end
            A
        end
    end
end


coefficients(f::AbstractProductFun)=funlist2coefficients(f.coefficients)

function coefficients{S,V,T<:Number}(f::AbstractProductFun{S,V,T},ox::FunctionSpace,oy::FunctionSpace)
    m=size(f,1)
    A=[pad!(coefficients(fx,ox),m) for fx in f.coefficients]
    B=hcat(A...)::Array{T,2}
    for k=1:size(B,1)
        B[k,:]=spaceconversion(vec(B[k,:]),space(f,2),oy)
    end
    
    B
end

# We assume that the spaces have the same values
function values(f::AbstractProductFun)
    n=size(f,1)
    M=hcat(map(fk->values(pad(fk,n)),f.coefficients)...)
    for k=1:size(f,1)
        M[k,:]=values(Fun(vec(M[k,:]),space(f,2)))
    end
    M     
end




points(f::AbstractProductFun,k)=points(space(f.coefficients[1]),size(f,k))

space(f::ProductFun,k::Integer)=k==1?error("ProductFun only has space(f,2)"):f.spacey
space(f::TensorFun,k::Integer)=k==1?space(f.coefficients[1]):f.spacey

domain(f::AbstractProductFun,k::Integer)=domain(space(f,k))



evaluate{S,V,T}(f::AbstractProductFun{S,V,T},x::Real,::Colon)=Fun(T[fc[x] for fc in f.coefficients],space(f,2))
evaluate(f::AbstractProductFun,x::Real,y::Real)=evaluate(f,x,:)[y]
evaluate(f::TensorFun,x::Colon,y::Real)=evaluate(f.',y,:)  # doesn't make sense For general product fon without specifying space
evaluate(f::AbstractProductFun,xx::Vector,yy::Vector)=hcat([evaluate(f,x,:)[[yy]] for x in xx]...).'
evaluate(f::AbstractProductFun,x::Range,y::Range)=evaluate(f,[x],[y])


*{F<:AbstractProductFun}(c::Number,f::F)=F(c*f.coefficients,f.spacey)
*(f::AbstractProductFun,c::Number)=c*f



chop{F<:AbstractProductFun}(f::F,es...)=F(map(g->chop(g,es...),f.coefficients),f.spacey)


##TODO: following assumes f is never changed....maybe should be deepcopy?
function +{F<:AbstractProductFun}(f::F,c::Number)
    cfs=copy(f.coefficients)
    cfs[1]+=c
    F(cfs,f.spacey)
end
+(c::Number,f::AbstractProductFun)=f+c
-(f::AbstractProductFun,c::Number)=f+(-c)
-(c::Number,f::AbstractProductFun)=c+(-f)


function +{F<:AbstractProductFun}(f::F,g::F)
    if size(f,2) >= size(g,2)
        @assert f.spacey==g.spacey
        cfs = copy(f.coefficients)
        for k=1:size(g,2)
            cfs[k]+=g.coefficients[k]
        end
        
        F(cfs,f.spacey)
    else
        g+f
    end
end

-(f::AbstractProductFun)=(-1)*f
-(f::AbstractProductFun,g::AbstractProductFun)=f+(-g)


LowRankFun(f::TensorFun)=LowRankFun(f.coefficients,space(f,2))

function differentiate(f::TensorFun,j::Integer)
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
    @eval ($op){S<:IntervalDomainSpace,V<:IntervalDomainSpace}(f::TensorFun{S,V})=TensorFun(chebyshevtransform($op(values(f))),domain(f))
end


for op = (:(Base.real),:(Base.imag),:(Base.conj)) 
    @eval ($op){S,V<:RealDomainSpace}(f::ProductFun{S,V}) = ProductFun(map($op,f.coefficients),f.spacey)
    @eval ($op){S,V<:RealDomainSpace}(f::TensorFun{S,V}) = TensorFun(map($op,f.coefficients),f.spacey)    
end

#For complex bases
Base.real{S,V,T}(u::TensorFun{S,V,T})=real(TensorFun(real(u.coefficients),space(u,2)).').'-imag(TensorFun(imag(u.coefficients),space(u,2)).').'
Base.imag{S,V,T}(u::TensorFun{S,V,T})=real(TensorFun(imag(u.coefficients),space(u,2)).').'+imag(TensorFun(real(u.coefficients),space(u,2)).').'








## ProductFun transform

function transform{ST<:FunctionSpace}(S::Vector{ST},T::FunctionSpace,V::Array)
    # We assume all S spaces have same domain/points
    C=Array(Complex{Float64},size(V)...)
    for k=1:size(V,1)
        C[k,:]=transform(T,vec(V[k,:]))
    end
    for k=1:size(C,2)
        C[:,k]=transform(S[k],C[:,k])
    end
    C
end

