
export ∂


abstract MultivariateDomain{T} <: Domain{T}
##TODO: MultivariateDomain{2}
abstract BivariateDomain{T} <: MultivariateDomain{T}





immutable ProductDomain{D<:Domain,T} <:BivariateDomain{T}
    domains::Vector{D}
end

ProductDomain{D<:Domain}(d::Vector{D})=ProductDomain{D,mapreduce(eltype,promote_type,d)}(d)

# product domains are their own canonical domain
for OP in (:fromcanonical,:tocanonical)
    @eval $OP(::ProductDomain,x,y)=(x,y)
end


ProductDomain(A,B)=ProductDomain([A,B])
*(A::Domain,B::Domain)=ProductDomain(A,B)

Base.length(d::ProductDomain)=length(d.domains)
Base.transpose(d::ProductDomain)=ProductDomain(d[2],d[1])
Base.getindex(d::ProductDomain,k::Integer)=d.domains[k]
==(d1::ProductDomain,d2::ProductDomain)=d1.domains==d2.domains

Base.first(d::ProductDomain)=(first(d[1]),first(d[2]))

function checkpoints(d::ProductDomain)
    ptsx=checkpoints(d[1])
    ptsy=checkpoints(d[2])
    ret=Array((Float64,Float64),0)
    for x in ptsx,y in ptsy
        push!(ret,(x,y))
    end
    ret
end

abstract MultivariateSpace{T} <: FunctionSpace{T}
abstract BivariateSpace{T} <: MultivariateSpace{T}


fromcanonical(d::MultivariateDomain,x::Tuple)=fromcanonical(d,x...)
tocanonical(d::MultivariateDomain,x::Tuple)=tocanonical(d,x...)
fromcanonical(d::MultivariateSpace,x...)=fromcanonical(domain(d),x...)
tocanonical(d::MultivariateSpace,x...)=tocanonical(domain(d),x...)


# This means x are represented as space S and y are represented as space T
abstract AbstractProductSpace{S,V,T} <: BivariateSpace{T}



immutable TensorSpace{S,V,T} <:AbstractProductSpace{S,V,T}
    spaces::(S,V)
end

for OP in (:spacescompatible,:(==))
    @eval $OP(A::TensorSpace,B::TensorSpace)=$OP(A.spaces[1],B.spaces[1])&&$OP(A.spaces[2],B.spaces[2])
end


TensorSpace(sp::Tuple)=TensorSpace{typeof(sp[1]),typeof(sp[2]),
                                                promote_type(eltype(sp[1]),eltype(sp[2])),
                                                Domain}(sp)


coefficient_type(S::TensorSpace,T)=promote_type(coefficient_type(S.spaces[1],T),coefficient_type(S.spaces[2],T))

TensorSpace(A,B)=TensorSpace((A,B))
TensorSpace(A::ProductDomain)=TensorSpace(Space(A[1]),Space(A[2]))
⊗(A::FunctionSpace,B::FunctionSpace)=TensorSpace(A,B)
domain(f::TensorSpace)=domain(f.spaces[1])*domain(f.spaces[2])
Space(sp::ProductDomain)=TensorSpace(sp)

*(A::FunctionSpace,B::FunctionSpace)=A⊗B


# every column is in the same space for a TensorSpace
columnspace(S::TensorSpace,::)=S.spaces[1]

Base.length(d::TensorSpace)=length(d.spaces)
Base.getindex(d::TensorSpace,k::Integer)=d.spaces[k]
=={S,V}(a::TensorSpace{S,V},b::TensorSpace{S,V})=(a.spaces[1]==b.spaces[1])&&(a.spaces[2]==b.spaces[2])

immutable ProductSpace{S<:FunctionSpace,T<:FunctionSpace} <: AbstractProductSpace{S,T}
    spacesx::Vector{S}
    spacey::T
end

coefficient_type(S::ProductSpace,T)=promote_type(coefficient_type(S.spacesx[1],T),coefficient_type(S.spacesy,T))

⊗{S<:FunctionSpace}(A::Vector{S},B::FunctionSpace)=ProductSpace(A,B)
domain(f::ProductSpace)=domain(f.spacesx[1])*domain(f.spacesy)

Base.getindex(d::ProductSpace,k::Integer)=k==1?d.spacesx:d.spacey


space(d::AbstractProductSpace,k)=d[k]
isambiguous(A::TensorSpace)=isambiguous(A[1])||isambiguous(A[2])

for TT in (:ProductDomain,:TensorSpace)
    @eval Base.transpose(d::$TT)=$TT(d[2],d[1])
end





##Transforms

function itransform!(S::TensorSpace,M::Matrix)
    n=size(M,1)
    for k=1:size(M,2)
        M[:,k]=itransform(space(S,1),M[:,k])
    end

    for k=1:n
        M[k,:]=itransform(space(S,2),vec(M[k,:]))
    end
    M
end

function itransform!(S::AbstractProductSpace,M::Matrix)
    n=size(M,1)

    ## The order matters
    pln=plan_itransform(columnspace(S,1),n)
    for k=1:size(M,2)
        M[:,k]=itransform(columnspace(S,k),M[:,k],pln)
    end

    for k=1:n
        M[k,:]=itransform(space(S,2),vec(M[k,:]))
    end
    M
end

function transform!(S::TensorSpace,M::Matrix)
    n=size(M,1)

    for k=1:size(M,2)
        M[:,k]=transform(space(S,1),M[:,k])
    end

    for k=1:n
        M[k,:]=transform(space(S,2),vec(M[k,:]))
    end
    M
end

function transform!{T}(S::AbstractProductSpace,M::Matrix{T})
    n=size(M,1)

    ## The order matters!!
    # For Disk Space, this is due to requiring decay
    # in function
    for k=1:n
        M[k,:]=transform(space(S,2),vec(M[k,:]))
    end

    pln=plan_transform(columnspace(S,1),n)
    for k=1:size(M,2)
        # col may not be full length
        col=transform(columnspace(S,k),M[:,k],pln)
        M[1:length(col),k]=col
        for j=length(col)+1:n
            M[j,k]=zero(T) # fill rest with zeros
        end
    end


    M
end



## points

points(d::Union(BivariateDomain,BivariateSpace),n,m)=points(d,n,m,1),points(d,n,m,2)

function points(d::BivariateSpace,n,m,k)
    ptsx=points(columnspace(d,1),n)
    ptst=points(space(d,2),m)

    promote_type(eltype(ptsx),eltype(ptst))[fromcanonical(d,x,t)[k] for x in ptsx, t in ptst]
end




##  Fun routines


function fromtensorind(k,j)
    n=k+j-2
    div(n*(n+1),2)+k
end

# which block of the tensor
# equivalent to sum of indices -1
totensorblock(n)=floor(Integer,sqrt(2n) + 1/2)
#gives the range corresponding to the block
fromtensorblock(j)=div(j*(j-1),2)+(1:j)

function totensorind(n)
    m=totensorblock(n)
    p=fromtensorind(1,m)
    j=1+n-p
    j,m-j+1
end


function fromtensor{T}(M::Matrix{T})
    ret=zeros(T,fromtensorind(size(M,1),size(M,2)))

    for k=1:size(M,1),j=1:size(M,2)
        ret[fromtensorind(k,j)]=M[k,j]
    end
    ret
end

function totensor{T}(M::Vector{T})
    inds=totensorind(length(M))
    m=inds[1]+inds[2]-1
    ret=zeros(T,m,m)
    for k=1:length(M)
        ret[totensorind(k)...]=M[k]
    end
    ret
end


function totree(v::Vector)
   m=totensorblock(length(v))
    r=Array(Vector{eltype(v)},m)
    for k=1:m-1
        r[k]=v[fromtensorblock(k)]
    end
    r[m]=pad!(v[fromtensorblock(m)[1]:end],m)
    r
end

fromtree{T}(v::Vector{Vector{T}})=vcat(v...)

function points(sp::TensorSpace,n)
    pts=Array((Float64,Float64),0)
    for x in points(sp[1],round(Int,sqrt(n))), y in points(sp[2],round(Int,sqrt(n)))
        push!(pts,(x,y))
    end
    pts
end

function transform!(S::TensorSpace,M::Matrix)
    n=size(M,1)

    for k=1:size(M,2)
        M[:,k]=transform(space(S,1),M[:,k])
    end

    for k=1:n
        M[k,:]=transform(space(S,2),vec(M[k,:]))
    end
    M
end

function transform(sp::TensorSpace,vals)
    m=round(Int,sqrt(length(vals)))
    M=reshape(vals,m,m)

    fromtensor(transform!(sp,M))
end

evaluate{S<:TensorSpace}(f::Fun{S},x)=ProductFun(totensor(coefficients(f)),space(f))[x...]
evaluate{S<:TensorSpace}(f::Fun{S},x,y)=ProductFun(totensor(coefficients(f)),space(f))[x,y]



coefficientmatrix{S<:AbstractProductSpace}(f::Fun{S})=totensor(f.coefficients)

Fun{T<:Number}(v::Vector{Vector{T}},S::TensorSpace)=Fun(fromtree(v),S)

## boundary


function ∂(d::ProductDomain)
    @assert length(d.domains) ==2
    #TODO: Generalize
    if isa(d[1],Interval)&&isa(d[2],Interval)
        PiecewiseInterval(d[1].a+im*d[2].a,d[1].b+im*d[2].a,d[1].b+im*d[2].b,d[1].a+im*d[2].b,d[1].a+im*d[2].a)
    elseif isa(d[1],Interval)&&isa(d[2],PeriodicInterval)
        UnionDomain([d[1].b+im*d[2],d[1].a+im*reverse(d[2])])
    elseif isa(d[1],PeriodicInterval)&&isa(d[2],Interval)
        UnionDomain([d[1]+im*d[2].a,reverse(d[1])+im*d[2].b])
    else
#        warn("∂ not implemented for "*string(typeof(d))*".  Returning [].")
        []
    end
end

#TODO: Implement
# function ∂(d::TensorSpace{Interval{Float64}})
#     @assert length(d.spaces) ==2
#     PiecewiseSpace([d[1].a+im*d[2],d[1].b+im*d[2],d[1]+im*d[2].a,d[1]+im*d[2].b])
# end
