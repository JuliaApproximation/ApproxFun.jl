

#Ultraspherical Spaces

type UltrasphericalSpace{T<:Union(IntervalDomain,DataType)} <: OperatorSpace
    order::Int
    domain::T     
end

UltrasphericalSpace(o::Integer)=UltrasphericalSpace(o,Any)

==(a::UltrasphericalSpace,b::UltrasphericalSpace)=a.order==b.order && a.domain==b.domain


##Check domain compatibility

domainscompatible(a::UltrasphericalSpace,b::UltrasphericalSpace) = a.domain == Any || b.domain == Any || a.domain == b.domain


spacescompatible(a::UltrasphericalSpace,b::UltrasphericalSpace) = domainscompatible(a,b) && a.order >= b.order

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


##TODO: Make general
function promoterangespace(P::Operator,sp::Space)
    psp = rangespace(P)
    
    if psp == Any || sp == Any || psp == sp
        P
    else
        @assert typeof(psp) <: UltrasphericalSpace
        @assert typeof(sp) <: UltrasphericalSpace
        @assert spacescompatible(sp,psp)
        
        if psp.order == sp.order
            P
        else
            ConversionOperator(psp.order:sp.order)*P
        end
    end
end




function promotedomainspace(P::Operator,sp::Space)
    psp = domainspace(P)
    
    if psp == Any || sp == Any || psp == sp
        P
    else
        @assert typeof(psp) <: UltrasphericalSpace
        @assert typeof(sp) <: UltrasphericalSpace
        @assert spacescompatible(psp,sp)
        
        
        if psp.order == sp.order
            P
        else
            P*ConversionOperator(sp.order:psp.order)
        end
    end
end


##TODO: remove?
promotespaces(op::BandedOperator,od::Range1)=promoterangespace(promotedomainspace(op,UltrasphericalSpace(od[1])),UltrasphericalSpace(od[end]))



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
