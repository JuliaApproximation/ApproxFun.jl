

#Ultraspherical Spaces

type UltrasphericalSpace{T<:Union(IntervalDomain,DataType)} <: OperatorSpace
    order::Int
    domain::T     
end

UltrasphericalSpace(o::Integer)=UltrasphericalSpace(o,Any)




##Check domain compatibility

domainscompatible(a::UltrasphericalSpace,b::UltrasphericalSpace) = a.domain == Any || b.domain == Any || a.domain == b.domain


#TODO: bad override?
==(a::UltrasphericalSpace,b::UltrasphericalSpace)=a.order==b.order && domainscompatible(a,b)

##max space



function Base.max(a::UltrasphericalSpace,b::UltrasphericalSpace)
    @assert domainscompatible(a,b)
    
    a.order > b.order?a:b
end

function Base.min(a::UltrasphericalSpace,b::UltrasphericalSpace)
    @assert domainscompatible(a,b)
    
    a.order < b.order?a:b
end



## Operator space manipulation



# DirichletSpaces


type UltrasphericalDirchletSpace{T<:Union(IntervalDomain,DataType)} <: OperatorSpace
    order::Int
    domain::T 
    left::Int
    right::Int    
end

==(a::UltrasphericalDirchletSpace,b::UltrasphericalDirchletSpace)=a.order==b.order && a.domain==b.domain && a.left==b.left && a.right==b.right


##Check domain compatibility

# domainscompatible(a::UltrasphericalDirchletSpace,b::UltrasphericalDirchletSpace) = a.domain == Any || b.domain == Any || a.domain == b.domain
# 
# 
# spacescompatible(a::UltrasphericalDirchletSpace,b::UltrasphericalSpace) = domainscompatible(a,b) && a.order >= b.order
