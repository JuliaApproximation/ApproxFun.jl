export Conversion


immutable Conversion{S<:FunctionSpace,V<:FunctionSpace,T<:Number} <: BandedOperator{T}
    domainspace::S
    rangespace::V
end
Conversion{S<:PeriodicDomainSpace,V<:PeriodicDomainSpace}(A::S,B::V)=Conversion{S,V,Complex{Float64}}(A,B)
Conversion{S<:IntervalDomainSpace,V<:IntervalDomainSpace}(A::S,B::V)=Conversion{S,V,Float64}(A,B)

domainspace(C::Conversion)=C.domainspace
rangespace(C::Conversion)=C.rangespace


Conversion(A::FunctionSpace,B::FunctionSpace,C::FunctionSpace)=Conversion(B,C)*Conversion(A,B)
function Conversion(a::FunctionSpace,b::FunctionSpace)
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
        Conversion{S,V}(a,b)
    end
end
    


## convert TO canonical
Conversion(A::FunctionSpace)=Conversion(A,canonicalspace(A))

