##
# ProductFun represents f(x,y) by Fun(.coefficients[k](x),.space[2])(y)
# where all coefficients are in the same space
##

export ProductFun

struct ProductFun{S<:UnivariateSpace,V<:UnivariateSpace,SS<:AbstractProductSpace,T}<:BivariateFun{T}
    coefficients::Vector{VFun{S,T}}     # coefficients are in x
    space::SS
end

ProductFun{S<:UnivariateSpace,V<:UnivariateSpace,T<:Number,DD<:BivariateDomain}(cfs::Vector{VFun{S,T}},sp::AbstractProductSpace{Tuple{S,V},DD}) =
    ProductFun{S,V,typeof(sp),T}(cfs,sp)
ProductFun{S<:UnivariateSpace,V<:UnivariateSpace,
         W<:UnivariateSpace,T<:Number,DD<:BivariateDomain}(cfs::Vector{VFun{S,T}},sp::AbstractProductSpace{Tuple{W,V},DD}) =
   ProductFun{W,V,typeof(sp),T}(VFun{W,T}[Fun(cfs[k],columnspace(sp,k)) for k=1:length(cfs)],sp)

Base.size(f::ProductFun,k::Integer) =
    k==1?mapreduce(ncoefficients,max,f.coefficients):length(f.coefficients)
Base.size(f::ProductFun) = (size(f,1),size(f,2))

## Construction in an AbstractProductSpace via a Matrix of coefficients

function ProductFun{S<:UnivariateSpace,V<:UnivariateSpace,T<:Number,DD<:BivariateDomain}(cfs::AbstractMatrix{T},sp::AbstractProductSpace{Tuple{S,V},DD};
                                                                       tol::Real=100eps(T),chopping::Bool=false)
    if chopping
        ncfs,kend=norm(cfs,Inf),size(cfs,2)
        if kend > 1 while isempty(chop(cfs[:,kend],ncfs*tol)) kend-=1 end end
        ret=VFun{S,T}[Fun(columnspace(sp,k),chop(cfs[:,k],ncfs*tol)) for k=1:max(kend,1)]
        ProductFun{S,V,typeof(sp),T}(ret,sp)
    else
        ret=VFun{S,T}[Fun(columnspace(sp,k),cfs[:,k]) for k=1:size(cfs,2)]
        ProductFun{S,V,typeof(sp),T}(ret,sp)
    end
end

## Construction in a ProductSpace via a Vector of Funs

function ProductFun{S<:UnivariateSpace,V<:UnivariateSpace,T<:Number}(M::Vector{VFun{S,T}},dy::V)
    funs=VFun{S,T}[Mk for Mk in M]
    ProductFun{S,V,ProductSpace{S,V},T}(funs,ProductSpace(S[space(fun) for fun in funs],dy))
end

## Adaptive construction

function ProductFun{S<:UnivariateSpace,V<:UnivariateSpace}(f::F,sp::AbstractProductSpace{Tuple{S,V}};tol=100eps())
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

function ProductFun(f::F,S::AbstractProductSpace,M::Integer,N::Integer;tol=100eps())
    xy = checkpoints(S)
    T = promote_type(eltype(f(first(xy)...)),eltype(S))
    ptsx,ptsy=points(S,M,N)
    vals=T[f(ptsx[k,j],ptsy[k,j]) for k=1:size(ptsx,1), j=1:size(ptsx,2)]
    ProductFun(transform!(S,vals),S;tol=tol,chopping=true)
end
ProductFun(f::F,S::TensorSpace) = ProductFun(LowRankFun(f,S))

ProductFun(f,dx::Space,dy::Space)=ProductFun(f,TensorSpace(dx,dy))

ProductFun(f::Function,dx::Space,dy::Space)=ProductFun(F(f),TensorSpace(dx,dy))
ProductFun(f::Function,args...;kwds...) = ProductFun(F(f),args...;kwds...)

## Domains promoted to Spaces

ProductFun(f::F,D::BivariateDomain,M::Integer,N::Integer)=ProductFun(f,Space(D),M,N)
ProductFun(f::F,d::Domain)=ProductFun(f,Space(d))
ProductFun(f::F,dx::UnivariateDomain,dy::UnivariateDomain)=ProductFun(f,Space(dx),Space(dy))
ProductFun(f::F) = ProductFun(f,Interval(),Interval())

## Conversion from other 2D Funs

ProductFun(f::LowRankFun)=ProductFun(coefficients(f),space(f,1),space(f,2))
ProductFun{S<:AbstractProductSpace}(f::Fun{S})=ProductFun(coefficientmatrix(f),space(f))

## Conversion to other ProductSpaces with the same coefficients

ProductFun(f::ProductFun,sp::TensorSpace)=space(f)==sp?f:ProductFun(coefficients(f,sp),sp)
ProductFun{S,V,SS<:TensorSpace}(f::ProductFun{S,V,SS},sp::ProductDomain)=ProductFun(f,Space(sp))

function ProductFun(f::ProductFun,sp::AbstractProductSpace)
    u=Array{VFun{typeof(columnspace(sp,1)),eltype(f)}}(length(f.coefficients))

    for k=1:length(f.coefficients)
        u[k]=Fun(f.coefficients[k],columnspace(sp,k))
    end

    ProductFun(u,sp)
end

## For specifying spaces by anonymous function

ProductFun(f::F,SF::Function,T::Space,M::Integer,N::Integer)=ProductFun(f,typeof(SF(1))[SF(k) for k=1:N],T,M)

## Conversion of a constant to a ProductFun

ProductFun(c::Number,sp::BivariateSpace)=ProductFun([Fun(c,columnspace(sp,1))],sp)
ProductFun(f::Fun,sp::BivariateSpace)=ProductFun([Fun(f,columnspace(sp,1))],sp)



## Utilities



function funlist2coefficients{S,T}(f::Vector{VFun{S,T}})
    A=zeros(T,mapreduce(ncoefficients,max,f),length(f))
    for k=1:length(f)
        A[1:ncoefficients(f[k]),k]=f[k].coefficients
    end
    A
end


function pad{S,V,SS,T}(f::ProductFun{S,V,SS,T},n::Integer,m::Integer)
    ret=Array{VFun{S,T}}(m)
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

function coefficients(f::ProductFun,ox::Space,oy::Space)
    T=eltype(f)
    m=size(f,1)
    B=Matrix{T}(m,length(f.coefficients))
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

(f::ProductFun)(x,y) = evaluate(f,x,y)

coefficients(f::ProductFun,ox::TensorSpace) = coefficients(f,ox[1],ox[2])




values{S,V,SS,T}(f::ProductFun{S,V,SS,T}) = itransform!(space(f),coefficients(f))


vecpoints{S,V,SS<:TensorSpace}(f::ProductFun{S,V,SS},k) = points(f.space[k],size(f,k))

space(f::ProductFun) = f.space
space(f::ProductFun,k) = space(space(f),k)
columnspace(f::ProductFun,k) = columnspace(space(f),k)

domain(f::ProductFun) = domain(f.space)
#domain(f::ProductFun,k)=domain(f.space,k)
canonicaldomain(f::ProductFun) = canonicaldomain(space(f))



function canonicalevaluate{S,V,SS,T}(f::ProductFun{S,V,SS,T},x::Number,::Colon)
    cd = canonicaldomain(f)
    Fun(setdomain(space(f,2),cd[2]),T[setdomain(fc,cd[1])(x) for fc in f.coefficients])
end
canonicalevaluate(f::ProductFun,x::Number,y::Number) = canonicalevaluate(f,x,:)(y)
canonicalevaluate{S,V,SS<:TensorSpace}(f::ProductFun{S,V,SS},x::Colon,y::Number) =
    evaluate(f.',y,:)  # doesn't make sense For general product fon without specifying space

canonicalevaluate(f::ProductFun,xx::AbstractVector,yy::AbstractVector) =
    hcat([evaluate(f,x,:)(yy) for x in xx]...).'


evaluate(f::ProductFun,x,y) = canonicalevaluate(f,tocanonical(f,x,y)...)

# TensorSpace does not use map
evaluate{S<:UnivariateSpace,V<:UnivariateSpace,SS<:TensorSpace,T}(f::ProductFun{S,V,SS,T},x::Number,::Colon) =
    Fun(factor(space(f),2),T[g(x) for g in f.coefficients])

evaluate{S<:UnivariateSpace,V<:UnivariateSpace,SS<:TensorSpace,T}(f::ProductFun{S,V,SS,T},x::Number,y::Number) =
    evaluate(f,x,:)(y)


evaluate(f::ProductFun,x) = evaluate(f,x...)

*{F<:ProductFun}(c::Number,f::F) = F(c*f.coefficients,f.space)
*(f::ProductFun,c::Number) = c*f


function chop{S}(f::ProductFun{S},es...)
    kend=size(f,2)
    if kend > 1 while isempty(chop(f.coefficients[kend].coefficients,es...)) kend-=1 end end
    ret=VFun{S,eltype(f)}[Fun(space(f.coefficients[k]),chop(f.coefficients[k].coefficients,es...)) for k=1:max(kend,1)]

    typeof(f)(ret,f.space)
end


##TODO: following assumes f is never changed....maybe should be deepcopy?
function +{F<:ProductFun}(f::F,c::Number)
    cfs=copy(f.coefficients)
    cfs[1]+=c
    F(cfs,f.space)
end
+(c::Number,f::ProductFun) = f+c
-(f::ProductFun,c::Number) = f+(-c)
-(c::Number,f::ProductFun) = c+(-f)


function +(f::ProductFun,g::ProductFun)
    if f.space == g.space
        if size(f,2) >= size(g,2)
            @assert f.space==g.space
            cfs = copy(f.coefficients)
            for k=1:size(g,2)
                cfs[k]+=g.coefficients[k]
            end

            ProductFun(cfs,f.space)
        else
            g+f
        end
    else
        s=conversion_type(f.space,g.space)
        ProductFun(f,s)+ProductFun(g,s)
    end
end

-(f::ProductFun) = (-1)*f
-(f::ProductFun,g::ProductFun) = f+(-g)

*(B::Fun,f::ProductFun) = ProductFun(map(c->B*c,f.coefficients),space(f))
*(f::ProductFun,B::Fun) = (B*f.').'


LowRankFun(f::ProductFun{S,V,SS}) where {S,V,SS<:TensorSpace} = LowRankFun(f.coefficients,factor(space(f),2))
LowRankFun(f::Fun) = LowRankFun(ProductFun(f))

function differentiate{S,V,SS<:TensorSpace}(f::ProductFun{S,V,SS},j::Integer)
    if j==1
        df=map(differentiate,f.coefficients)
        ProductFun(df,space(first(df)),factor(space(f),2))
    else
        differentiate(f.',1).'
    end
end

# If the transpose of the space exists, then the transpose of the ProductFun exists
Base.transpose{S,V,SS,T}(f::ProductFun{S,V,SS,T}) =
    ProductFun(transpose(coefficients(f)),transpose(space(f)))





for op in (:(Base.sin),:(Base.cos))
    @eval ($op)(f::ProductFun) =
        Fun(space(f),transform!(space(f),$op(values(pad(f,size(f,1)+20,size(f,2))))))
end

^(f::ProductFun,k::Integer) =
    Fun(space(f),transform!(space(f),values(pad(f,size(f,1)+20,size(f,2))).^k))

for op = (:(Base.real),:(Base.imag),:(Base.conj))
    @eval ($op)(f::ProductFun{S,V,SS}) where {S,V<:RealSpace,SS<:TensorSpace} =
        ProductFun(map($op,f.coefficients),space(f))
end

#For complex bases
Base.real{S,V,SS<:TensorSpace}(f::ProductFun{S,V,SS}) =
    real(ProductFun(real(u.coefficients),space(u)).').'-imag(ProductFun(imag(u.coefficients),space(u)).').'
Base.imag{S,V,SS<:TensorSpace}(f::ProductFun{S,V,SS}) =
    real(ProductFun(imag(u.coefficients),space(u)).').'+imag(ProductFun(real(u.coefficients),space(u)).').'



## Call LowRankFun version
# TODO: should cumsum and integrate return TensorFun or lowrankfun?
for op in (:(Base.sum),:(Base.cumsum),:integrate)
    @eval $op{S,V,SS<:TensorSpace}(f::ProductFun{S,V,SS},n...) = $op(LowRankFun(f),n...)
end


## ProductFun transform

# function transform{ST<:Space,N<:Number}(::Type{N},S::Vector{ST},T::Space,V::AbstractMatrix)
#     @assert length(S)==size(V,2)
#     # We assume all S spaces have same domain/points
#     C=Vector{N}(size(V)...)
#     for k=1:size(V,1)
#         C[k,:]=transform(T,vec(V[k,:]))
#     end
#     for k=1:size(C,2)
#         C[:,k]=transform(S[k],C[:,k])
#     end
#     C
# end
# transform{ST<:Space,N<:Real}(S::Vector{ST},T::Space{Float64},V::AbstractMatrix{N})=transform(Float64,S,T,V)
# transform{ST<:Space}(S::Vector{ST},T::Space,V::AbstractMatrix)=transform(Complex{Float64},S,T,V)




for op in (:tocanonical,:fromcanonical)
    @eval $op(f::ProductFun,x...) = $op(space(f),x...)
end
