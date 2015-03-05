## Drop space drops the first n entries from a space

immutable SliceSpace{index,stride,DS,T,D}<: FunctionSpace{T,D}
    space::DS

    SliceSpace(sp::DS)=new(sp)
    SliceSpace(d::Domain)=new(DS(d))
end


spacescompatible{n,st,DS,T,D}(S1::SliceSpace{n,st,DS,T,D},S2::SliceSpace{n,st,DS,T,D})=spacescompatible(S1.space,S2.space)

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


getindex{ind,DS,T,D}(E::Evaluation{SliceSpace{ind,1,DS,T,D},Bool},kr::Range)=Evaluation(E.space.space,E.x,E.order)[kr+ind]
getindex{ind,DS,T,D}(E::Evaluation{SliceSpace{ind,1,DS,T,D}},kr::Range)=Evaluation(E.space.space,E.x,E.order)[kr+ind]

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
        a  # we can write droping coefficients as a banded operator
    else
        NoSpace()
    end
end



## Resolve conflict
for TYP in (:ReImSpace,:ReSpace,:ImSpace)
    @eval begin
        coefficients(::Vector,sp::$TYP,slp::SliceSpace)=error("coefficients not implemented from "*string(typeof(sp))*" to "*string(typeof(slp)))
        coefficients(::Vector,sp::SliceSpace,slp::$TYP)=error("coefficients not implemented from "*typeof(sp)*" to "*typeof(slp))
    end
end

# v[k]=v[stride*k+index]
function coefficients(v::Vector,sp::SliceSpace,dropsp::SliceSpace)
    if sp==dropsp
        v
    else
        coefficients(v,sp,canonicalspace(sp),dropsp)
    end
end

function coefficients(v::Vector,sp::FunctionSpace,dropsp::SliceSpace)
    if sp==dropsp.space
        n=index(dropsp)
        st=stride(dropsp)
        v[st+n:st:end]
    else
        coefficients(v,sp,canonicalspace(dropsp),dropsp)
    end
end

function coefficients{V}(v::Vector{V},dropsp::SliceSpace,sp::FunctionSpace)
    if sp==dropsp.space
        n=index(dropsp)
        st=stride(dropsp)
        ret=zeros(V,st*length(v)+n)
        ret[st+n:st:end]=v
        ret
    else
        coefficients(v,dropsp,canonicalspace(dropsp),sp)
    end
end

canonicalspace(a::SliceSpace)=a.space

## points


## transform
# TODO: padding shouldn't be necessary
function transform{n,st,S,T}(sp::SliceSpace{n,st,S},vals::Vector{T})
    ret=transform(sp.space,vals)
    # TODO: Test for zeros?
    [ret[st+n:st:end],zeros(T,n)] # ensure same length by padding with zeros
end

itransform{n,st,S,T}(sp::SliceSpace{n,st,S},cfs::Vector{T})=itransform(sp.space,coefficients(cfs,sp,sp.space))

##TODO: coefficients

points{n,S,T}(f::Fun{SliceSpace{n,1,S,T}})=points(space(f),length(f)+n)

## ProductFUn

values{S<:SliceSpace}(f::ProductFun{S})=values(ProductFun(f,space(f,1).space,space(f,2)))

function coefficients{n,DS,TT,D}(f::ProductFun{SliceSpace{n,1,DS,TT,D}},ox::FunctionSpace,oy::FunctionSpace)
    T=eltype(f)
    m=size(f,1)
    A=[pad!(coefficients(fx,ox),m+n) for fx in f.coefficients]
    B=hcat(A...)::Array{T,2}
    for k=1:size(B,1)
        ccfs=coefficients(vec(B[k,:]),space(f,2),oy)
        if length(ccfs)>size(B,2)
            B=pad(B,size(B,1),length(ccfs))
        end
        B[k,1:length(ccfs)]=ccfs
        #B[k,length(ccfs):1:end]=zero(T)
    end

    B
end
