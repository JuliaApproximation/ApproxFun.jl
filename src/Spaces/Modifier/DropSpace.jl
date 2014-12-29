## Drop space drops the first n entries from a space

immutable StrideSpace{index,stride,DS,T,D}<: FunctionSpace{T,D}
    space::DS 
    
    StrideSpace(sp::DS)=new(sp)
    StrideSpace(d::Domain)=new(DS(d))
end

#typealias DropSpace{n,DS,T,D} StrideSpace{n,1,DS,T,D}

StrideSpace{T,D}(sp::FunctionSpace{T,D},n::Integer,st::Integer)=StrideSpace{n,st,typeof(sp),T,D}(sp)
StrideSpace{T,D}(sp::FunctionSpace{T,D},n::Integer)=StrideSpace(sp,n,1)

domain(DS::StrideSpace)=domain(DS.space)
bandinds{n,st,S,T,D}(C::Conversion{StrideSpace{n,st,S,T,D},S})=-n,0

function addentries!{ind,st,S,T,D}(C::Conversion{StrideSpace{ind,st,S,T,D},S},A,kr::Range)
    ds =domainspace(C)
    @assert st==1

    for k=max(kr[1],ind+1):kr[end]
        A[k,k-ind]+=1
    end
    A
end

=={n,st,S,T,D}(a::StrideSpace{n,st,S,T,D},b::StrideSpace{n,st,S,T,D})=a.space==b.space

function conversion_rule{n,S<:FunctionSpace,T,D}(a::StrideSpace{n,1,S,T,D},b::StrideSpace{n,1,S,T,D})
     @assert a==b
     a
end
# return the space that has banded Conversion to the other
function conversion_rule{n,S<:FunctionSpace,T,D}(a::StrideSpace{n,1,S,T,D},b::S)
    @assert a.space==b
    a
end



## Resolve conflict
spaceconversion{n,st,S1<:FunctionSpace,S2<:FunctionSpace,T1,T2,V,D}(f::Vector{V},a::ReImSpace{S1,T1},b::StrideSpace{n,st,S2,T2,D})=error("Not implemented")
spaceconversion{n,st,S1<:FunctionSpace,S2<:FunctionSpace,T1,T2,V,D}(f::Vector{V},b::StrideSpace{n,st,S2,T2,D},a::ReImSpace{S1,T1})=error("Not implemented")
# v[k]=v[stride*k+index]
spaceconversion{n,st,S1<:FunctionSpace,S2<:StrideSpace,U,V,D}(v::Vector{V},sp::StrideSpace{n,st,S1,U,D},dropsp::StrideSpace{n,st,S2,U,D})=error("Not implemented")
function spaceconversion{n,st,S1<:FunctionSpace,S2<:FunctionSpace,U,V,D}(v::Vector{V},sp::StrideSpace{n,st,S1,U,D},dropsp::StrideSpace{n,st,S2,U,D})
    if sp==dropsp
        v
    else
        error("spaceconversion not implemented from "*typeof(sp))
    end
end

function spaceconversion{n,st,S<:FunctionSpace,U,V,D}(v::Vector{V},sp::S,dropsp::StrideSpace{n,st,S,U,D})
    @assert sp==dropsp.space
    v[st+n:st:end]
end

function spaceconversion{n,st,S<:FunctionSpace,U,V,D}(v::Vector{V},dropsp::StrideSpace{n,st,S,U,D},sp::S)
    @assert sp==dropsp.space

    ret=zeros(V,st*length(v)+n)
    ret[st+n:st:end]=v
    ret
end

canonicalspace(a::StrideSpace)=a.space




## transform
# TODO: padding shouldn't be necessary
function transform{n,st,S,T}(sp::StrideSpace{n,st,S},vals::Vector{T})
    ret=transform(sp.space,vals)
    # TODO: Test for zeros?
    [ret[st+n:st:end],zeros(T,n)] # ensure same length by padding with zeros
end

itransform{n,st,S,T}(sp::StrideSpace{n,st,S},cfs::Vector{T})=itransform(sp.space,spaceconversion(cfs,sp,sp.space))

##TODO: spaceconversion



