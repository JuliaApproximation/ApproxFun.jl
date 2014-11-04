## Drop space drops the first n entries from a space

immutable DropSpace{S,n,T}<: DomainSpace{T}
    space::S 
end

DropSpace{T}(sp::DomainSpace{T},n::Integer)=DropSpace{typeof(sp),n,T}(sp)

domain(DS::DropSpace)=domain(DS.space)
bandinds{S,n,T}(::Conversion{DropSpace{S,n,T},S})=-n,0

function addentries!{S,T,n}(C::Conversion{DropSpace{S,n,T},S},A::ShiftArray,kr::Range)
    for k=max(kr[1],n+1):kr[end]
        A[k,-n]+=1
    end
    A
end

=={S,n,T}(a::DropSpace{S,n,T},b::DropSpace{S,n,T})=a.space==b.space

function conversion_rule{S<:FunctionSpace,n,T}(a::DropSpace{S,n,T},b::DropSpace{S,n,T})
    @assert a==b
    a
end
# return the space that has banded Conversion to the other
function conversion_rule{S<:FunctionSpace,n,T}(a::DropSpace{S,n,T},b::S)
    @assert a.space==b
    a
end


canonicalspace(a::DropSpace)=a.space



## transform

function transform{S,n}(sp::DropSpace{S,n},vals)
    tol=1e-12

    ret=transform(sp.space,vals)
    @assert norm(ret[1:n])<tol
    
    ret[n+1:end]
end

itransform{S,n,T}(sp::DropSpace{S,n},cfs::Vector{T})=itransform(sp.space,[zeros(T,n),cfs])