## Drop space drops the first n entries from a space

immutable SliceSpace{index,stride,DS,T,D}<: FunctionSpace{T,D}
    space::DS 
    
    SliceSpace(sp::DS)=new(sp)
    SliceSpace(d::Domain)=new(DS(d))
end


index{n}(::SliceSpace{n})=n
Base.stride{n,st}(::SliceSpace{n,st})=st

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
     if a==b
        a
    else
        NoSpace()
    end
end
# return the space that has banded Conversion to the other
function conversion_rule{n,S<:FunctionSpace,T,D}(a::SliceSpace{n,1,S,T,D},b::FunctionSpace)
    if a.space==b
        b  # we can write droping coefficients as a banded operator
    else
        NoSpace()
    end
end



## Resolve conflict
spaceconversion(::Vector,sp::ReImSpace,slp::SliceSpace)=error("spaceconversion not implemented from "*typeof(sp)*" to "*typeof(slp))
spaceconversion(::Vector,sp::SliceSpace,slp::ReImSpace)=error("spaceconversion not implemented from "*typeof(sp)*" to "*typeof(slp))

# v[k]=v[stride*k+index]
function spaceconversion(v::Vector,sp::SliceSpace,dropsp::SliceSpace)
    if sp==dropsp
        v
    else
        error("spaceconversion not implemented from "*typeof(sp))
    end
end

function spaceconversion(v::Vector,sp::FunctionSpace,dropsp::SliceSpace)
    @assert sp==dropsp.space
    n=index(dropsp)
    st=stride(dropsp)
    v[st+n:st:end]
end

function spaceconversion{V}(v::Vector{V},dropsp::SliceSpace,sp::FunctionSpace)
    @assert sp==dropsp.space
    n=index(dropsp)
    st=stride(dropsp)
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



