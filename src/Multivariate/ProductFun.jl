##
# ProductFun represents f(x,y) by Fun(.coefficients[k][x],.space[2])[y]
# where all coefficients are in the same space
##

export ProductFun

immutable ProductFun{S<:UnivariateSpace,V<:UnivariateSpace,SS<:AbstractProductSpace,T}<:BivariateFun
    coefficients::Vector{Fun{S,T}}     # coefficients are in x
    space::SS
end

ProductFun{S<:UnivariateSpace,V<:UnivariateSpace,T<:Number,P}(cfs::Vector{Fun{S,T}},sp::AbstractProductSpace{@compat(Tuple{S,V}),P,2})=ProductFun{S,V,typeof(sp),T}(cfs,sp)
ProductFun{S<:UnivariateSpace,V<:UnivariateSpace,T<:Number,P}(cfs::Vector{Fun{S,T}},sp::TensorSpace{@compat(Tuple{S,V}),P,2})=ProductFun{S,V,typeof(sp),T}(cfs,sp)  #for ambiguity

function ProductFun{S<:UnivariateSpace,V<:UnivariateSpace,
                    W<:UnivariateSpace,T<:Number,P}(cfs::Vector{Fun{S,T}},sp::AbstractProductSpace{@compat(Tuple{W,V}),P,2})
   ProductFun{W,V,typeof(sp),T}(Fun{W,T}[Fun(cfs[k],columnspace(sp,k)) for k=1:length(cfs)],sp)
end

#for ambiguity
function ProductFun{S<:UnivariateSpace,V<:UnivariateSpace,
                    W<:UnivariateSpace,T<:Number,P}(cfs::Vector{Fun{S,T}},sp::TensorSpace{@compat(Tuple{W,V}),P,2})
   ProductFun{W,V,typeof(sp),T}(Fun{W,T}[Fun(cfs[k],columnspace(sp,k)) for k=1:length(cfs)],sp)
end

Base.size(f::ProductFun,k::Integer)=k==1?mapreduce(length,max,f.coefficients):length(f.coefficients)
Base.size(f::ProductFun)=(size(f,1),size(f,2))
Base.eltype{S,V,SS,T}(::ProductFun{S,V,SS,T})=T

## Construction in an AbstractProductSpace via a Matrix of coefficients

function ProductFun{S<:UnivariateSpace,V<:UnivariateSpace,T<:Number,P}(cfs::Matrix{T},sp::AbstractProductSpace{@compat(Tuple{S,V}),P,2};
                                                                        tol::Real=100eps(T),chopping::Bool=false)
    if chopping
        ncfs,kend=norm(cfs,Inf),size(cfs,2)
        if kend > 1 while isempty(chop(cfs[:,kend],ncfs*tol)) kend-=1 end end
        ret=Fun{S,T}[Fun(chop(cfs[:,k],ncfs*tol),columnspace(sp,k)) for k=1:max(kend,1)]
        ProductFun{S,V,typeof(sp),T}(ret,sp)
    else
        ret=Fun{S,T}[Fun(cfs[:,k],columnspace(sp,k)) for k=1:size(cfs,2)]
        ProductFun{S,V,typeof(sp),T}(ret,sp)
    end
end


#repeated for ambiguity
function ProductFun{S<:UnivariateSpace,V<:UnivariateSpace,T<:Number,P}(cfs::Matrix{T},sp::TensorSpace{@compat(Tuple{S,V}),P,2};
                                                                        tol::Real=100eps(T),chopping::Bool=false)
    if chopping
        ncfs,kend=norm(cfs,Inf),size(cfs,2)
        if kend > 1 while isempty(chop(cfs[:,kend],ncfs*tol)) kend-=1 end end
        ret=Fun{S,T}[Fun(chop(cfs[:,k],ncfs*tol),columnspace(sp,k)) for k=1:max(kend,1)]
        ProductFun{S,V,typeof(sp),T}(ret,sp)
    else
        ret=Fun{S,T}[Fun(cfs[:,k],columnspace(sp,k)) for k=1:size(cfs,2)]
        ProductFun{S,V,typeof(sp),T}(ret,sp)
    end
end

## Construction in a ProductSpace via a Vector of Funs

function ProductFun{S<:UnivariateSpace,V<:UnivariateSpace,T<:Number}(M::Vector{Fun{S,T}},dy::V)
    funs=Fun{S,T}[Mk for Mk in M]
    ProductFun{S,V,ProductSpace{S,V},T}(funs,ProductSpace(S[space(fun) for fun in funs],dy))
end

## Adaptive construction

function ProductFun{S<:UnivariateSpace,V<:UnivariateSpace,T}(f,sp::AbstractProductSpace{@compat(Tuple{S,V}),T,2};tol=100eps())
    for n = 50:100:5000
        X = coefficients(ProductFun(f,sp,n,n;tol=tol))
        if size(X,1)<n && size(X,2)<n
            return ProductFun(X,sp;tol=tol)
        end
    end
    warn("Maximum grid size of ("*string(5000)*","*string(5000)*") reached")
    ProductFun(f,sp,5000,5000;tol=tol,chopping=true)
end

## ProductFun values to coefficients

function ProductFun(f,S::AbstractProductSpace,M::Integer,N::Integer;tol=100eps())
    xy = checkpoints(S)
    T = promote_type(eltype(f(first(xy)...)),eltype(S))
    ptsx,ptsy=points(S,M,N)
    vals=T[f(ptsx[k,j],ptsy[k,j]) for k=1:size(ptsx,1), j=1:size(ptsx,2)]
    ProductFun(transform!(S,vals),S;tol=tol,chopping=true)
end
ProductFun{SS<:UnivariateSpace,V<:UnivariateSpace,P}(f,S::TensorSpace{@compat(Tuple{SS,V}),P,2}) = ProductFun(LowRankFun(f,S))
ProductFun(f,dx::FunctionSpace,dy::FunctionSpace)=ProductFun(f,TensorSpace(dx,dy))


## Domains promoted to FunctionSpaces

ProductFun(f,D::BivariateDomain,M::Integer,N::Integer)=ProductFun(f,Space(D),M,N)
ProductFun(f,d::Domain)=ProductFun(f,Space(d))
ProductFun(f,dx::UnivariateDomain,dy::UnivariateDomain)=ProductFun(f,Space(dx),Space(dy))
ProductFun(f) = ProductFun(f,Interval(),Interval())

## Conversion from other 2D Funs

ProductFun(f::LowRankFun)=ProductFun(coefficients(f),space(f,1),space(f,2))
ProductFun{S<:AbstractProductSpace}(f::Fun{S})=ProductFun(coefficientmatrix(f),space(f))

## Conversion to other ProductSpaces with the same coefficients

ProductFun{S<:UnivariateSpace,V<:UnivariateSpace,P}(f::ProductFun,sp::TensorSpace{@compat(Tuple{S,V}),P,2})=space(f)==sp?f:ProductFun(coefficients(f,sp),sp)
ProductFun{S<:UnivariateSpace,V<:UnivariateSpace,P}(f::ProductFun,sp::AbstractProductSpace{@compat(Tuple{S,V}),P,2})=space(f)==sp?f:ProductFun(coefficients(f,sp),sp)
ProductFun{S,V,SS<:TensorSpace}(f::ProductFun{S,V,SS},sp::ProductDomain)=ProductFun(f,Space(sp))

## For specifying spaces by anonymous function

ProductFun(f,SF,T::FunctionSpace,M::Integer,N::Integer)=ProductFun(f,typeof(SF(1))[SF(k) for k=1:N],T,M)

## Conversion of a constant to a ProductFun

ProductFun{S<:UnivariateSpace,V<:UnivariateSpace,P}(c::Union(Fun,Number),sp::TensorSpace{@compat(Tuple{S,V}),P,2})=ProductFun([Fun(c,columnspace(sp,1))],sp)  # for ambiguity
ProductFun{S<:UnivariateSpace,V<:UnivariateSpace,P}(c::Union(Fun,Number),sp::AbstractProductSpace{@compat(Tuple{S,V}),P,2})=ProductFun([Fun(c,columnspace(sp,1))],sp)  # for ambiguity
ProductFun(c::Union(Fun,Number),sp::BivariateSpace)=ProductFun([Fun(c,columnspace(sp,1))],sp)




## Utilities



function funlist2coefficients{S,T}(f::Vector{Fun{S,T}})
    A=zeros(T,mapreduce(length,max,f),length(f))
    for k=1:length(f)
        A[1:length(f[k]),k]=f[k].coefficients
    end
    A
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

function pad!{S,V,SS,T}(f::ProductFun{S,V,SS,T},::Colon,m::Integer)
    cm=length(f.coefficients)
    resize!(f.coefficients,m)

    for k=cm+1:m
        f.coefficients[k]=zero(T,columnspace(f,k))
    end
    f
end


coefficients(f::ProductFun)=funlist2coefficients(f.coefficients)

function coefficients(f::ProductFun,ox::FunctionSpace,oy::FunctionSpace)
    T=eltype(f)
    m=size(f,1)
    B=Array(T,m,length(f.coefficients))
    # convert in x direction
    #TODO: adaptively grow in x?
    for k=1:length(f.coefficients)
        B[:,k]=pad!(coefficients(f.coefficients[k],ox),m)
    end
    
    # convert in y direction
    for k=1:size(B,1)
        ccfs=coefficients(vec(B[k,:]),space(f,2),oy)
        if length(ccfs)>size(B,2)
            B=pad(B,size(B,1),length(ccfs))
        end
        B[k,1:length(ccfs)]=ccfs
        for j=length(ccfs)+1:size(B,2)
            B[k,j]=zero(T)
        end
    end

    B
end

coefficients(f::ProductFun,ox::TensorSpace)=coefficients(f,ox[1],ox[2])

values{S,V,SS,T}(f::ProductFun{S,V,SS,T})=itransform!(space(f),coefficients(f))


vecpoints{S,V,SS<:TensorSpace}(f::ProductFun{S,V,SS},k)=points(f.space[k],size(f,k))

space(f::ProductFun)=f.space
space(f::ProductFun,k)=space(space(f),k)
columnspace(f::ProductFun,k)=columnspace(space(f),k)

domain(f::ProductFun)=domain(f.space)
#domain(f::ProductFun,k)=domain(f.space,k)




canonicalevaluate{S,V,SS,T}(f::ProductFun{S,V,SS,T},x::Number,::Colon)=Fun(T[fc[x] for fc in f.coefficients],space(f,2))
canonicalevaluate(f::ProductFun,x::Number,y::Number)=canonicalevaluate(f,x,:)[y]
canonicalevaluate{S,V,SS<:TensorSpace}(f::ProductFun{S,V,SS},x::Colon,y::Number)=evaluate(f.',y,:)  # doesn't make sense For general product fon without specifying space

canonicalevaluate(f::ProductFun,xx::Vector,yy::Vector)=hcat([evaluate(f,x,:)[[yy]] for x in xx]...).'


evaluate(f::ProductFun,x,y)=canonicalevaluate(f,tocanonical(f,x,y)...)
evaluate(f::ProductFun,x::Range,y::Range)=evaluate(f,[x],[y])


*{F<:ProductFun}(c::Number,f::F)=F(c*f.coefficients,f.space)
*(f::ProductFun,c::Number)=c*f


#.'
function chop{S}(f::ProductFun{S},es...)
    kend=size(f,2)
    if kend > 1 while isempty(chop(f.coefficients[kend].coefficients,es...)) kend-=1 end end
    ret=Fun{S,eltype(f)}[Fun(chop(f.coefficients[k].coefficients,es...),space(f.coefficients[k])) for k=1:max(kend,1)]

    typeof(f)(ret,f.space)
end


##TODO: following assumes f is never changed....maybe should be deepcopy?
function +{F<:ProductFun}(f::F,c::Number)
    cfs=copy(f.coefficients)
    cfs[1]+=c
    F(cfs,f.space)
end
+(c::Number,f::ProductFun)=f+c
-(f::ProductFun,c::Number)=f+(-c)
-(c::Number,f::ProductFun)=c+(-f)


function +{F<:ProductFun}(f::F,g::F)
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

-(f::ProductFun)=(-1)*f
-(f::ProductFun,g::ProductFun)=f+(-g)

*(B::Fun,f::ProductFun)=ProductFun(map(c->B*c,f.coefficients),space(f))
*(f::ProductFun,B::Fun)=(B*f.').'


LowRankFun{S,V,SS<:TensorSpace}(f::ProductFun{S,V,SS})=LowRankFun(f.coefficients,space(f,2))

function differentiate{S,V,SS<:TensorSpace}(f::ProductFun{S,V,SS},j::Integer)
    if j==1
        df=map(differentiate,f.coefficients)
        ProductFun(df,space(first(df)),space(f,2))
    else
        differentiate(f.',1).'
    end
end

# If the transpose of the space exists, then the transpose of the ProductFun exists
Base.transpose{S,V,SS,T}(f::ProductFun{S,V,SS,T})=ProductFun(transpose(coefficients(f)),transpose(space(f)))





for op in (:(Base.sin),:(Base.cos))
    @eval ($op)(f::ProductFun)=Fun(transform!(space(f),$op(values(pad(f,size(f,1)+20,size(f,2))))),space(f))
end

.^(f::ProductFun,k::Integer)=Fun(transform!(space(f),values(pad(f,size(f,1)+20,size(f,2))).^k),space(f))

for op = (:(Base.real),:(Base.imag),:(Base.conj))
    @eval ($op){S,V<:FunctionSpace{RealBasis},SS<:TensorSpace}(f::ProductFun{S,V,SS}) = ProductFun(map($op,f.coefficients),space(f))
end

#For complex bases
Base.real{S,V,SS<:TensorSpace}(f::ProductFun{S,V,SS})=real(ProductFun(real(u.coefficients),space(u)).').'-imag(ProductFun(imag(u.coefficients),space(u)).').'
Base.imag{S,V,SS<:TensorSpace}(f::ProductFun{S,V,SS})=real(ProductFun(imag(u.coefficients),space(u)).').'+imag(ProductFun(real(u.coefficients),space(u)).').'



## Call LowRankFun version
# TODO: should cumsum and integrate return TensorFun or lowrankfun?
for op in (:(Base.sum),:(Base.cumsum),:integrate)
    @eval $op{S,V,SS<:TensorSpace}(f::ProductFun{S,V,SS},n...)=$op(LowRankFun(f),n...)
end


## ProductFun transform

# function transform{ST<:FunctionSpace,N<:Number}(::Type{N},S::Vector{ST},T::FunctionSpace,V::Matrix)
#     @assert length(S)==size(V,2)
#     # We assume all S spaces have same domain/points
#     C=Array(N,size(V)...)
#     for k=1:size(V,1)
#         C[k,:]=transform(T,vec(V[k,:]))
#     end
#     for k=1:size(C,2)
#         C[:,k]=transform(S[k],C[:,k])
#     end
#     C
# end
# transform{ST<:FunctionSpace,N<:Real}(S::Vector{ST},T::FunctionSpace{Float64},V::Matrix{N})=transform(Float64,S,T,V)
# transform{ST<:FunctionSpace}(S::Vector{ST},T::FunctionSpace,V::Matrix)=transform(Complex{Float64},S,T,V)




for op in (:tocanonical,:fromcanonical)
    @eval $op(f::ProductFun,x...)=$op(space(f),x...)
end



