export ProductFun



##
# ProductFun represents f(x,y) by Fun(.coefficients[k][x],.space[2])[y]
# where all coefficients are in the same space
##


immutable ProductFun{S<:FunctionSpace,V<:FunctionSpace,SS<:AbstractProductSpace,T}<:BivariateFun
    coefficients::Vector{Fun{S,T}}     # coefficients are in x
    space::SS
end

typealias TensorFun{S,V,T} ProductFun{S,V,TensorSpace{S,V},T}


ProductFun{S<:FunctionSpace,V<:FunctionSpace,T<:Number}(cfs::Vector{Fun{S,T}},sp::AbstractProductSpace{S,V})=ProductFun{S,V,typeof(sp),T}(cfs,sp)




function ProductFun{S,T}(M::Vector{Fun{S,T}},dy::FunctionSpace)
    Sx=typeof(M[1].space)
    funs=Fun{Sx,T}[Mk for Mk in M]
    ProductFun{Sx,typeof(dy),ProductSpace{Sx,typeof(dy)},$T}(funs,ProductSpace(Sx[space(fun) for fun in funs],dy))
end

ProductFun(M,dx::FunctionSpace,dy::FunctionSpace)=ProductFun(M,TensorSpace(dx,dy))

function ProductFun{T<:Number,S<:FunctionSpace,V<:FunctionSpace}(cfs::Matrix{T},D::AbstractProductSpace{S,V})
     ret=Array(Fun{S,T},size(cfs,2))
     for k=1:size(cfs,2)
         ret[k]=chop!(Fun(cfs[:,k],columnspace(D,k)),10eps())
     end
     ProductFun{S,V,typeof(D),T}(ret,D)
end



ProductFun(f::Function,dy::Domain)=error("This function is only implemented to avoid ambiguity, do not call.")
ProductFun(f::Function,d::BivariateDomain)=ProductFun(f,Space(d))
ProductFun(f,dy::Domain)=ProductFun(f,Space(dy))
ProductFun(f,d::BivariateDomain)=ProductFun(f,Space(d))
ProductFun(f,dx::Domain,dy::Domain)=ProductFun(f,Space(dx),Space(dy))
ProductFun(f::LowRankFun)=ProductFun(coefficients(f),space(f,1),space(f,2))

#Need to templt TensorFun because parameters are in subfields
ProductFun(f::ProductFun,sp::AbstractProductSpace)=space(f)==sp?f:ProductFun(coefficients(f,sp),sp)

ProductFun{S,V}(f::TensorFun{S,V},sp1::Domain,sp2::Domain)=ProductFun(f,Space(sp1),Space(sp2))
ProductFun{S,V}(f::TensorFun{S,V},sp::ProductDomain)=ProductFun(f,Space(sp))

ProductFun(f::Function,d1...)=ProductFun(LowRankFun(f,d1...))
# function ProductFun{ST<:FunctionSpace}(f::Function,S::Vector{ST},T::FunctionSpace,N::Integer)
#     M=length(S)     # We assume the list of spaces is the same length
#     xx=points(S[1],N)
#     tt=points(T,M)
#     V=Float64[f(x,t) for x=xx, t=tt] #TODO:type
#     ProductFun(transform(S,T,V),S,T)
# end

function ProductFun{SS<:FunctionSpace{Float64},VV<:FunctionSpace{Float64}}(f::Function,S::AbstractProductSpace{SS,VV},N::Integer,M::Integer)
    ptsx,ptsy=points(S,N,M)
    V=Float64[f(ptsx[k,j],ptsy[k,j]) for k=1:size(ptsx,1), j=1:size(ptsx,2)]#TODO:type
    ProductFun(transform!(S,V),S)
end

function ProductFun(f::Function,S::AbstractProductSpace,N::Integer,M::Integer)
    ptsx,ptsy=points(S,N,M)

    f1=f(ptsx[1,1],ptsy[1,1])
    T=coefficient_type(S,typeof(f1))  # type of coefficients
    # use coefficient type so we can reuse array
    V=T[f(ptsx[k,j],ptsy[k,j]) for k=1:size(ptsx,1), j=1:size(ptsx,2)]
    ProductFun(transform!(S,V),S)
end

ProductFun(f::Function,D::BivariateDomain,N::Integer,M::Integer)=ProductFun(f,Space(D),N,M)

ProductFun(f::Function,D::TensorSpace)=ProductFun(LowRankFun(f,D))
ProductFun(f::Function,D::ProductDomain)=ProductFun(LowRankFun(f,D))

function ProductFun(f::Function,D)
     Nmax=400

     tol=1E-12

     for N=50:25:Nmax
         X=coefficients(ProductFun(f,D,N,N))
         if norm(X[end-3:end,:])<tol && norm(X[:,end-3:end])<tol
             chop!(X,tol)
             return ProductFun(X,D)
         end
     end
     error("Maximum grid reached")
     ProductFun(f,D,Nmax,Nmax)
end





# For specifying spaces by anonymous function
ProductFun(f::Function,SF::Function,T::FunctionSpace,N::Integer,M::Integer)=ProductFun(f,typeof(SF(1))[SF(k) for k=1:M],T,N)


# convert a constant  to a 2D Fun

ProductFun(c::Number,sp::BivariateFunctionSpace)=ProductFun([Fun(c,columnspace(sp,1))],sp)
ProductFun(f::Fun,sp::BivariateFunctionSpace)=ProductFun([Fun(f,columnspace(sp,1))],sp)



Base.size(f::ProductFun,k::Integer)=k==1?mapreduce(length,max,f.coefficients):length(f.coefficients)
Base.size(f::ProductFun)=(size(f,1),size(f,2))
Base.eltype{S,V,SS,T}(::ProductFun{S,V,SS,T})=T

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
    A=[pad!(coefficients(fx,ox),m) for fx in f.coefficients]
    B=hcat(A...)::Array{T,2}
    for k=1:size(B,1)
        ccfs=spaceconversion(vec(B[k,:]),space(f,2),oy)
        if length(ccfs)>size(B,2)
            B=pad(B,size(B,1),length(ccfs))
        end
        B[k,1:length(ccfs)]=ccfs
        B[k,length(ccfs):1:end]=zero(T)
    end

    B
end

coefficients(f::ProductFun,ox::TensorSpace)=coefficients(f,ox[1],ox[2])

values{S,V,SS,T}(f::ProductFun{S,V,SS,T})=itransform!(space(f),coefficients(f))


vecpoints{S,V,T}(f::TensorFun{S,V,T},k)=points(f.space[k],size(f,k))

space(f::ProductFun)=f.space
space(f::ProductFun,k)=space(space(f),k)
columnspace(f::ProductFun,k)=columnspace(space(f),k)

domain(f::ProductFun)=domain(f.space)
#domain(f::ProductFun,k)=domain(f.space,k)



canonicalevaluate{S,V,SS,T}(f::ProductFun{S,V,SS,T},x::Real,::Colon)=Fun(T[fc[x] for fc in f.coefficients],space(f,2))
canonicalevaluate(f::ProductFun,x::Real,y::Real)=canonicalevaluate(f,x,:)[y]
canonicalevaluate{S,V}(f::TensorFun{S,V},x::Colon,y::Real)=evaluate(f.',y,:)  # doesn't make sense For general product fon without specifying space
canonicalevaluate(f::ProductFun,xx::Vector,yy::Vector)=hcat([evaluate(f,x,:)[[yy]] for x in xx]...).'


evaluate(f::ProductFun,x,y)=canonicalevaluate(f,tocanonical(f,x,y)...)
evaluate(f::ProductFun,x::Range,y::Range)=evaluate(f,[x],[y])


*{F<:ProductFun}(c::Number,f::F)=F(c*f.coefficients,f.space)
*(f::ProductFun,c::Number)=c*f



chop{F<:ProductFun}(f::F,es...)=F(map(g->chop(g,es...),f.coefficients),f.space)


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


LowRankFun{S,V}(f::TensorFun{S,V})=LowRankFun(f.coefficients,space(f,2))

function differentiate{S,V}(f::TensorFun{S,V},j::Integer)
    if j==1
        df=map(differentiate,f.coefficients)
        ProductFun(df,space(first(df)),space(f,2))
    else
        differentiate(f.',1).'
    end
end


Base.transpose{S,V}(f::TensorFun{S,V})=ProductFun(coefficients(f).',space(f,2),space(f,1))





for op in (:(Base.sin),:(Base.cos))
    @eval ($op)(f::ProductFun)=Fun(transform!(space(f),$op(values(pad(f,size(f,1)+20,size(f,2))))),space(f))
end

.^(f::ProductFun,k::Integer)=Fun(transform!(space(f),values(pad(f,size(f,1)+20,size(f,2))).^k),space(f))

for op = (:(Base.real),:(Base.imag),:(Base.conj))
    @eval ($op){S,V<:FunctionSpace{RealBasis}}(f::TensorFun{S,V}) = ProductFun(map($op,f.coefficients),space(f))
end

#For complex bases
Base.real{S,V,T}(u::TensorFun{S,V,T})=real(ProductFun(real(u.coefficients),space(u)).').'-imag(ProductFun(imag(u.coefficients),space(u)).').'
Base.imag{S,V,T}(u::TensorFun{S,V,T})=real(ProductFun(imag(u.coefficients),space(u)).').'+imag(ProductFun(real(u.coefficients),space(u)).').'



## Call LowRankFun verison
# TODO: should cumsum and integrate return TensorFun or lowrankfun?
for op in (:(Base.sum),:(Base.cumsum),:integrate)
    @eval $op{S,V}(f::TensorFun{S,V},n...)=$op(LowRankFun(f),n...)
end


# Let K be a ProductFun and f be a Fun.
# op(f,K) acts as operating in the x variable, and
# op(K,f) acts as operating in the y variable.

#=
for op = (:*,:.*,:./,:/)
    #@eval ($op){S,T,U,V}(f::Fun{S,T},A::Vector{Fun{U,V}})=map(a->($op)(f,a),A)
    @eval ($op)(f::Fun,K::ProductFun) = ProductFun(($op)(f,K.A),K.B)
    #@eval ($op){S,T,U,V}(B::Vector{Fun{U,V}},f::Fun{S,T})=map(b->($op)(b,f),B)
    @eval ($op)(K::ProductFun,f::Fun) = ProductFun(K.A,($op)(K.B,f))
end
=#

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

include("ChebyshevAddition.jl")
