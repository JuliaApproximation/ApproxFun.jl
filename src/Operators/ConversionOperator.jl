export ConversionOperator


immutable ConversionOperator{S<:FunctionSpace,V<:FunctionSpace} <: BandedOperator{Float64}
    domainspace::S
    rangespace::V
end

domainspace(C::ConversionOperator)=C.domainspace
rangespace(C::ConversionOperator)=C.rangespace


ConversionOperator(A::FunctionSpace,B::FunctionSpace,C::FunctionSpace)=ConversionOperator(B,C)*ConversionOperator(A,B)
function ConversionOperator(a::FunctionSpace,b::FunctionSpace)
    if a==b
        IdentityOperator()
    elseif conversion_type(a,b)==NoSpace()
        sp=canonicalspace(a)
        if typeof(sp) == typeof(a)
            error("implement ConversionOperator from " * string(typeof(sp)) * " to " * string(typeof(b)))
        elseif typeof(sp) == typeof(b)
            error("implement ConversionOperator from " * string(typeof(a)) * " to " * string(typeof(sp)))
        else        
            ConversionOperator(a,sp,b)
        end
    else
        ConversionOperator{S,V}(a,b)
    end
end
    


## convert TO canonical
ConversionOperator(A::FunctionSpace)=ConversionOperator(A,canonicalspace(A))

