## Drop space drops the first n entries from a space

immutable SliceSpace{index,stride,DS,T,DD,dim}<: Space{T,DD,dim}
    space::DS

    SliceSpace(sp::DS)=new(sp)
    SliceSpace(d::Domain{Number,dim})=new(DS(d))
end


spacescompatible{n,st,DS,T,DD,d}(S1::SliceSpace{n,st,DS,T,DD,d},S2::SliceSpace{n,st,DS,T,DD,d})=spacescompatible(S1.space,S2.space)

index{n}(::SliceSpace{n})=n
Base.stride{n,st}(::SliceSpace{n,st})=st

SliceSpace(sp::Space,n::Integer,st::Integer)=SliceSpace{n,st,typeof(sp),basistype(sp),domaintype(sp),ndims(sp)}(sp)
SliceSpace(sp,n::Integer)=SliceSpace(sp,n,1)

domain(DS::SliceSpace)=domain(DS.space)
bandinds{n,st,S,T,DD,d}(C::Conversion{SliceSpace{n,st,S,T,DD,d},S})=-n,0

function addentries!{ind,st,S,T,DD,d}(C::Conversion{SliceSpace{ind,st,S,T,DD,d},S},A,kr::Range,::Colon)
    ds =domainspace(C)
    @assert st==1

    for k=max(kr[1],ind+1):kr[end]
        A[k,k-ind]+=1
    end
    A
end


getindex{ind,DS,T,DD,d}(E::Evaluation{SliceSpace{ind,1,DS,T,DD,d},Bool},kr::Range)=Evaluation(E.space.space,E.x,E.order)[kr+ind]
getindex{ind,DS,T,DD,d}(E::Evaluation{SliceSpace{ind,1,DS,T,DD,d}},kr::Range)=Evaluation(E.space.space,E.x,E.order)[kr+ind]

=={n,st,S,T,DD,d}(a::SliceSpace{n,st,S,T,DD,d},b::SliceSpace{n,st,S,T,DD,d})=a.space==b.space

function conversion_rule{n,S<:Space,T}(a::SliceSpace{n,1,S,T},b::SliceSpace{n,1,S,T})
     if a==b
        a
    else
        NoSpace()
    end
end
# return the space that has banded Conversion to the other
function conversion_rule{n,S<:Space,T}(a::SliceSpace{n,1,S,T},b::Space)
    if a.space==b
        a  # we can write droping coefficients as a banded operator
    else
        NoSpace()
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

values{S<:SliceSpace,V<:SliceSpace}(f::ProductFun{S,V})=values(ProductFun(f,space(f,1).space,space(f,2).space))
values{S<:SliceSpace}(f::ProductFun{S})=values(ProductFun(f,space(f,1).space,space(f,2)))


function coefficients{n,DS,TT,DD}(f::ProductFun{SliceSpace{n,1,DS,TT,DD,1}},ox::Space,oy::Space)
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
