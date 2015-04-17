
export TensorSpace,ProductSpace

#  SV is a tuple of d spaces
abstract AbstractProductSpace{SV,T,d} <: FunctionSpace{T,d}


immutable TensorSpace{SV,T,d} <:AbstractProductSpace{SV,T,d}
    spaces::SV
end

for OP in (:spacescompatible,:(==))
    @eval $OP{SV,T,d}(A::TensorSpace{SV,T,d},B::TensorSpace{SV,T,d})=all(Bool[$OP(A.spaces[k],B.spaces[k]) for k=1:length(A.spaces)])
end


TensorSpace(sp::Tuple)=TensorSpace{typeof(sp),mapreduce(basistype,promote_type,sp),mapreduce(ndims,+,sp)}(sp)


coefficient_type(S::TensorSpace,T)=mapreduce(sp->coefficient_type(sp,T),promote_type,S.spaces)

TensorSpace(A...)=TensorSpace(tuple(A...))
TensorSpace(A::ProductDomain)=TensorSpace(tuple(map(Space,A.domains)...))
⊗{SV1,T1,SV2,T2,d1,d2}(A::TensorSpace{SV1,T1,d1},B::TensorSpace{SV2,T2,d2})=TensorSpace(A.spaces...,B.spaces...)
⊗{SV,T,V,d1,d2}(A::TensorSpace{SV,T,d1},B::FunctionSpace{V,d2})=TensorSpace(A.spaces...,B)
⊗{SV,T,V,d1,d2}(A::FunctionSpace{V,d2},B::TensorSpace{SV,T,d1})=TensorSpace(A,B.spaces...)
⊗{T,V,d1,d2}(A::FunctionSpace{T,d1},B::FunctionSpace{V,d2})=TensorSpace(A,B)

domain(f::TensorSpace)=mapreduce(domain,*,f.spaces)
Space(sp::ProductDomain)=TensorSpace(sp)

*(A::FunctionSpace,B::FunctionSpace)=A⊗B


# every column is in the same space for a TensorSpace
#TODO: remove
columnspace(S::TensorSpace,::)=S.spaces[1]

Base.length(d::TensorSpace)=length(d.spaces)
Base.getindex(d::TensorSpace,k::Integer)=d.spaces[k]


immutable ProductSpace{S<:FunctionSpace,V<:FunctionSpace,T} <: AbstractProductSpace{(S,V),T,2}
    spacesx::Vector{S}
    spacey::V
end

ProductSpace(spacesx::Vector,spacey)=ProductSpace{eltype(spacesx),
                                                  typeof(spacey),
                                                  promote_type(basistype(first(spacesx)),basistype(spacey))}(spacesx,spacey)

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


#TODO: Implement
# function ∂(d::TensorSpace{Interval{Float64}})
#     @assert length(d.spaces) ==2
#     PiecewiseSpace([d[1].a+im*d[2],d[1].b+im*d[2],d[1]+im*d[2].a,d[1]+im*d[2].b])
# end
