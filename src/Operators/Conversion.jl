export Conversion

abstract AbstractConversion{T}<:BandedOperator{T}

immutable Conversion{S<:FunctionSpace,V<:FunctionSpace,T<:Number} <: AbstractConversion{T}
    domainspace::S
    rangespace::V
end
# Conversion{S<:PeriodicDomainSpace,V<:PeriodicDomainSpace}(A::S,B::V)=Conversion{S,V,Complex{Float64}}(A,B)
# Conversion{S<:IntervalDomainSpace,V<:IntervalDomainSpace}(A::S,B::V)=Conversion{S,V,Float64}(A,B)

domainspace(C::Conversion)=C.domainspace
rangespace(C::Conversion)=C.rangespace


Conversion(A::FunctionSpace,B::FunctionSpace,C::FunctionSpace)=Conversion(B,C)*Conversion(A,B)

#TODO: Periodic
function Conversion{T,V}(a::DomainSpace{T},b::DomainSpace{V})
    if a==b
        IdentityOperator()
    elseif conversion_type(a,b)==NoSpace()
        sp=canonicalspace(a)
        if typeof(sp) == typeof(a)
            error("implement Conversion from " * string(typeof(sp)) * " to " * string(typeof(b)))
        elseif typeof(sp) == typeof(b)
            error("implement Conversion from " * string(typeof(a)) * " to " * string(typeof(sp)))
        else        
            Conversion(a,sp,b)
        end
    else
        Conversion{typeof(a),typeof(b),promote_type(T,V)}(a,b)
    end
end
    


## convert TO canonical
Conversion(A::FunctionSpace)=Conversion(A,canonicalspace(A))



## Wrapper
# this allows for a Derivative implementation to return another operator, use a SpaceOperator containing
# the domain and range space
# but continue to know its a derivative

type ConversionWrapper{S<:BandedOperator} <: AbstractConversion{Float64}
    op::S
end

addentries!(D::ConversionWrapper,A::ShiftArray,k::Range)=addentries!(D.op,A,k)
for func in (:rangespace,:domainspace,:bandinds)
    @eval $func(D::ConversionWrapper)=$func(D.op)
end
