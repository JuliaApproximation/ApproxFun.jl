export TensorFun,ProductFun


abstract AbstractProductFun{S,V,T} <: BivariateFun


immutable TensorFun{S<:FunctionSpace,V<:FunctionSpace,T<:Union(Float64,Complex{Float64})}<:AbstractProductFun{S,V,T}
    coefficients::Vector{Fun{S,T}}     # coefficients are in x
    space::TensorSpace{S,V}
end

immutable ProductFun{S<:FunctionSpace,V<:FunctionSpace,SS<:AbstractProductSpace{S,V},T<:Union(Float64,Complex{Float64})}<:AbstractProductFun{S,V,T}
    coefficients::Vector{Fun{S,T}}     # coefficients are in x
    space::SS
end

ProductFun{S<:FunctionSpace,V<:FunctionSpace,T<:Number}(cfs::Vector{Fun{S,T}},sp::AbstractProductSpace{S,V})=ProductFun{S,V,typeof(sp),T}(cfs,sp)


for T in (:Float64,:(Complex{Float64}))
    @eval begin
        TensorFun{S}(M::Vector{Fun{S,$T}},dy::FunctionSpace)=TensorFun{typeof(M[1].space),typeof(dy),$T}(Fun{typeof(M[1].space),$T}[Mk for Mk in M],TensorSpace(M[1].space,dy))
        function ProductFun{S}(M::Vector{Fun{S,$T}},dy::FunctionSpace)
            Sx=typeof(M[1].space)
            funs=Fun{Sx,$T}[Mk for Mk in M]
            ProductFun{Sx,typeof(dy),ProductSpace{Sx,typeof(dy)},$T}(funs,ProductSpace(Sx[space(fun) for fun in funs],dy))    
        end
    end
end


function TensorFun{T<:Number,S<:FunctionSpace}(cfs::Matrix{T},dx::S,dy::FunctionSpace)
    ret=Array(Fun{S,T},size(cfs,2))
    for k=1:size(cfs,2)
        ret[k]=chop!(Fun(cfs[:,k],dx),10eps())
    end
    TensorFun{S,typeof(dy),T}(ret,TensorSpace(dx,dy))
end

function ProductFun{T<:Number,S<:FunctionSpace,V<:FunctionSpace}(cfs::Matrix{T},D::AbstractProductSpace{S,V})
     ret=Array(Fun{S,T},size(cfs,2))
     for k=1:size(cfs,2)
         ret[k]=chop!(Fun(cfs[:,k],columnspace(D,k)),10eps())
     end
     ProductFun{S,V,typeof(D),T}(ret,D)
end

TensorFun{T<:Number}(cfs::Matrix{T},d::TensorSpace)=TensorFun(cfs,d[1],d[2])
function TensorFun{F<:Fun}(cfs::Vector{F},d::TensorSpace)
    for cf in cfs
        @assert space(cf)==d[1]
    end
    
    TensorFun(F,d[2])
end
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
    V=Float64[f(x,t) for x=xx, t=tt] #TODO:type
    ProductFun(transform(S,T,V),S,T)
end

function ProductFun(f::Function,S::AbstractProductSpace,N::Integer,M::Integer)
    ptsx,ptsy=points(S,N,M)
    V=Float64[f(ptsx[k,j],ptsy[k,j]) for k=1:size(ptsx,1), j=1:size(ptsx,2)]#TODO:type
    ProductFun(transform(S,V),S)
end

ProductFun(f::Function,D::BivariateDomain,N::Integer,M::Integer)=ProductFun(f,Space(D),N,M)



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


function pad{S,V,T}(f::TensorFun{S,V,T},n::Integer,m::Integer)
    ret=Array(Fun{S,T},m)
    cm=min(length(f.coefficients),m)
    for k=1:cm
        ret[k]=pad(f.coefficients[k],n)
    end
    zr=zero(space(f,1))
    for k=cm+1:m
        ret[k]=zr
    end
    TensorFun{S,V,T}(ret,f.space)
end

function pad{S,V,SS,T}(f::ProductFun{S,V,SS,T},n::Integer,m::Integer)
    ret=Array(Fun{S,T},m)
    cm=min(length(f.coefficients),m)
    for k=1:cm
        ret[k]=pad(f.coefficients[k],n)
    end

    for k=cm+1:m
        ret[k]=zero(T,columnspace(f,k))
    end
    ProductFun{S,V,SS,T}(ret,f.space)
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
function values{S,V,T}(f::TensorFun{S,V,T})
    M=coefficients(f)
    n=size(M,1)
    for k=1:size(M,2)
        M[:,k]=itransform(space(f,1),M[:,k])
    end

    n=size(M,1)
    for k=1:n
        M[k,:]=itransform(space(f,2),vec(M[k,:]))
    end 
    M  
end


function values(f::ProductFun)
    M=coefficients(f)
    n=size(M,1)
    pts=points(columnspace(f,1),n)
    for k=1:size(M,2)
        M[:,k]=itransform(columnspace(f,k),M[:,k],pts)
    end

    n=size(M,1)
    for k=1:n
        M[k,:]=itransform(space(f,2),vec(M[k,:]))
    end 
    M  
end




points(f::TensorFun,k)=points(space(f,k),size(f,k))
points(f::ProductFun,k...)=points(f.space,size(f,1),size(f,2),k...)


space(f::AbstractProductFun)=f.space
space(f::AbstractProductFun,k)=space(space(f),k)
columnspace(f::ProductFun,k)=columnspace(space(f),k)

domain(f::AbstractProductFun)=domain(f.space)
#domain(f::AbstractProductFun,k)=domain(f.space,k)



canonicalevaluate{S,V,T}(f::AbstractProductFun{S,V,T},x::Real,::Colon)=Fun(T[fc[x] for fc in f.coefficients],space(f,2))
canonicalevaluate(f::AbstractProductFun,x::Real,y::Real)=canonicalevaluate(f,x,:)[y]
canonicalevaluate(f::TensorFun,x::Colon,y::Real)=evaluate(f.',y,:)  # doesn't make sense For general product fon without specifying space
canonicalevaluate(f::AbstractProductFun,xx::Vector,yy::Vector)=hcat([evaluate(f,x,:)[[yy]] for x in xx]...).'


evaluate(f::ProductFun,x,y)=canonicalevaluate(f,tocanonical(f,x,y)...)
evaluate(f::TensorFun,x,y)=canonicalevaluate(f,x,y)


evaluate(f::TensorFun,x::Range,y::Range)=canonicalevaluate(f,[x],[y])
evaluate(f::ProductFun,x::Range,y::Range)=evaluate(f,[x],[y])


*{F<:AbstractProductFun}(c::Number,f::F)=F(c*f.coefficients,f.space)
*(f::AbstractProductFun,c::Number)=c*f



chop{F<:AbstractProductFun}(f::F,es...)=F(map(g->chop(g,es...),f.coefficients),f.space)


##TODO: following assumes f is never changed....maybe should be deepcopy?
function +{F<:AbstractProductFun}(f::F,c::Number)
    cfs=copy(f.coefficients)
    cfs[1]+=c
    F(cfs,f.space)
end
+(c::Number,f::AbstractProductFun)=f+c
-(f::AbstractProductFun,c::Number)=f+(-c)
-(c::Number,f::AbstractProductFun)=c+(-f)


function +{F<:AbstractProductFun}(f::F,g::F)
    if size(f,2) >= size(g,2)
        @assert f.space==g.space
        cfs = copy(f.coefficients)
        for k=1:size(g,2)
            cfs[k]+=g.coefficients[k]
        end
        
        F(cfs,f.space)
    else
        g+f
    end
end

-(f::AbstractProductFun)=(-1)*f
-(f::AbstractProductFun,g::AbstractProductFun)=f+(-g)


LowRankFun(f::TensorFun)=LowRankFun(f.coefficients,space(f,2))

function differentiate(f::TensorFun,j::Integer)
    if j==1
        TensorFun(map(diff,f.coefficients),space(f,2))
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
#    @eval ($op){S,V<:DomainSpace{Flaot64}}(f::ProductFun{S,V}) = ProductFun(map($op,f.coefficients),space(f,2))
    @eval ($op){S,V<:DomainSpace{Float64}}(f::TensorFun{S,V}) = TensorFun(map($op,f.coefficients),space(f,2))    
end

#For complex bases
Base.real{S,V,T}(u::TensorFun{S,V,T})=real(TensorFun(real(u.coefficients),space(u,2)).').'-imag(TensorFun(imag(u.coefficients),space(u,2)).').'
Base.imag{S,V,T}(u::TensorFun{S,V,T})=real(TensorFun(imag(u.coefficients),space(u,2)).').'+imag(TensorFun(real(u.coefficients),space(u,2)).').'








## ProductFun transform

function transform{ST<:FunctionSpace,N<:Number}(::Type{N},S::Vector{ST},T::FunctionSpace,V::Matrix)
    @assert length(S)==size(V,2)
    # We assume all S spaces have same domain/points
    C=Array(N,size(V)...)
    for k=1:size(V,1)
        C[k,:]=transform(T,vec(V[k,:]))
    end
    for k=1:size(C,2)
        C[:,k]=transform(S[k],C[:,k])
    end
    C
end
transform{ST<:FunctionSpace,N<:Real}(S::Vector{ST},T::DomainSpace{Float64},V::Matrix{N})=transform(Float64,S,T,V)
transform{ST<:FunctionSpace}(S::Vector{ST},T::FunctionSpace,V::Matrix)=transform(Complex{Float64},S,T,V)




for op in (:tocanonical,:fromcanonical)
    @eval $op(f::AbstractProductFun,x...)=$op(space(f),x...)
end

