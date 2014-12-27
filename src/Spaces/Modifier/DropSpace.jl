## Drop space drops the first n entries from a space

immutable DropSpace{DS,n,T,D}<: FunctionSpace{T,D}
    space::DS 
    DropSpace(sp::DS)=new(sp)
    DropSpace(d::Domain)=new(DS(d))
end

DropSpace{T,D}(sp::FunctionSpace{T,D},n::Integer)=DropSpace{typeof(sp),n,T,D}(sp)

domain(DS::DropSpace)=domain(DS.space)
bandinds{S,n,T,D}(::Conversion{DropSpace{S,n,T,D},S})=-n,0

function addentries!{S,T,n,D}(C::Conversion{DropSpace{S,n,T,D},S},A,kr::Range)
    for k=max(kr[1],n+1):kr[end]
        A[k,k-n]+=1
    end
    A
end

=={S,n,T,D}(a::DropSpace{S,n,T,D},b::DropSpace{S,n,T,D})=a.space==b.space

function conversion_rule{S<:FunctionSpace,n,T,D}(a::DropSpace{S,n,T,D},b::DropSpace{S,n,T,D})
    @assert a==b
    a
end
# return the space that has banded Conversion to the other
function conversion_rule{S<:FunctionSpace,n,T,D}(a::DropSpace{S,n,T,D},b::S)
    @assert a.space==b
    a
end



## Resolve conflict
function spaceconversion{S1<:FunctionSpace,S2<:FunctionSpace,n,T1,T2,V,D}(f::Vector{V},a::ReImSpace{S1,T1},b::DropSpace{S2,n,T2,D})
     error("Not implemented")
end


function spaceconversion{S<:FunctionSpace,n,U,V,D}(v::Vector{V},sp::S,dropsp::DropSpace{S,n,U,D})
    @assert sp==dropsp.space
    @assert norm(v[1:n])<100eps()
    v[n+1:end]
end

canonicalspace(a::DropSpace)=a.space




## transform
# TODO: padding shouldn't be necessary
function transform{S,n,T}(sp::DropSpace{S,n},vals::Vector{T})
    tol=1e-12

    ret=transform(sp.space,vals)
    @assert norm(ret[1:n])<tol
    
    [ret[n+1:end],zeros(T,n)] # ensure same length by padding with zeros
end

itransform{S,n,T}(sp::DropSpace{S,n},cfs::Vector{T})=itransform(sp.space,[zeros(T,n),cfs[1:end-n]])


##TODO: spaceconversion



