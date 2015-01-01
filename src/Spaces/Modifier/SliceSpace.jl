## Drop space drops the first n entries from a space

immutable SliceSpace{index,stride,DS,T,D}<: FunctionSpace{T,D}
    space::DS 
    
    SliceSpace(sp::DS)=new(sp)
    SliceSpace(d::Domain)=new(DS(d))
end

#typealias DropSpace{n,DS,T,D} SliceSpace{n,1,DS,T,D}

SliceSpace{T,D}(sp::FunctionSpace{T,D},n::Integer,st::Integer)=SliceSpace{n,st,typeof(sp),T,D}(sp)
SliceSpace{T,D}(sp::FunctionSpace{T,D},n::Integer)=SliceSpace(sp,n,1)

domain(DS::SliceSpace)=domain(DS.space)
bandinds{n,st,S,T,D}(C::Conversion{SliceSpace{n,st,S,T,D},S})=-n,0

function addentries!{ind,st,S,T,D}(C::Conversion{SliceSpace{ind,st,S,T,D},S},A,kr::Range)
    ds =domainspace(C)
    @assert st==1

    for k=max(kr[1],ind+1):kr[end]
        A[k,k-ind]+=1
    end
    A
end

=={n,st,S,T,D}(a::SliceSpace{n,st,S,T,D},b::SliceSpace{n,st,S,T,D})=a.space==b.space

function conversion_rule{n,S<:FunctionSpace,T,D}(a::SliceSpace{n,1,S,T,D},b::SliceSpace{n,1,S,T,D})
     @assert a==b
     a
end
# return the space that has banded Conversion to the other
function conversion_rule{n,S<:FunctionSpace,T,D}(a::SliceSpace{n,1,S,T,D},b::S)
    @assert a.space==b
    a
end



## Resolve conflict
spaceconversion{n,st,S1<:FunctionSpace,S2<:FunctionSpace,T1,T2,V,D}(f::Vector{V},a::ReImSpace{S1,T1},b::SliceSpace{n,st,S2,T2,D})=error("Not implemented")
spaceconversion{n,st,S1<:FunctionSpace,S2<:FunctionSpace,T1,T2,V,D}(f::Vector{V},b::SliceSpace{n,st,S2,T2,D},a::ReImSpace{S1,T1})=error("Not implemented")
# v[k]=v[stride*k+index]
spaceconversion{n,st,S1<:FunctionSpace,S2<:SliceSpace,U,V,D}(v::Vector{V},sp::SliceSpace{n,st,S1,U,D},dropsp::SliceSpace{n,st,S2,U,D})=error("Not implemented")
function spaceconversion{n,st,S1<:FunctionSpace,S2<:FunctionSpace,U,V,D}(v::Vector{V},sp::SliceSpace{n,st,S1,U,D},dropsp::SliceSpace{n,st,S2,U,D})
    if sp==dropsp
        v
    else
        error("spaceconversion not implemented from "*typeof(sp))
    end
end

function spaceconversion{n,st,S<:FunctionSpace,U,V,D}(v::Vector{V},sp::S,dropsp::SliceSpace{n,st,S,U,D})
    @assert sp==dropsp.space
    v[st+n:st:end]
end

function spaceconversion{n,st,S<:FunctionSpace,U,V,D}(v::Vector{V},dropsp::SliceSpace{n,st,S,U,D},sp::S)
    @assert sp==dropsp.space

    ret=zeros(V,st*length(v)+n)
    ret[st+n:st:end]=v
    ret
end

canonicalspace(a::SliceSpace)=a.space




## transform
# TODO: padding shouldn't be necessary
function transform{n,st,S,T}(sp::SliceSpace{n,st,S},vals::Vector{T})
    ret=transform(sp.space,vals)
    # TODO: Test for zeros?
    [ret[st+n:st:end],zeros(T,n)] # ensure same length by padding with zeros
end

itransform{n,st,S,T}(sp::SliceSpace{n,st,S},cfs::Vector{T})=itransform(sp.space,spaceconversion(cfs,sp,sp.space))

##TODO: spaceconversion



