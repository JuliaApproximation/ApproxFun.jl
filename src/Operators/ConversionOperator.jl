export ConversionOperator

ConversionOperator(::ChebyshevSpace,::ChebyshevSpace)=IdentityOperator() #TODO: should this be disallowed
ConversionOperator(A::FunctionSpace,B::FunctionSpace,C::FunctionSpace)=ConversionOperator(B,C)*ConversionOperator(A,B)
function ConversionOperator(a::FunctionSpace,b::FunctionSpace)
    if a==b
        IdentityOperator()
    elseif conversion_type(a,b) != NoSpace()
        error("Override ConversionOperator if you override conversion_type(" * string(typeof(a)) *","*string(typeof(b))*")") 
    elseif typeof(a) <: ChebyshevSpace
        error("Override ConversionOperator from ChebyshevSpace to " * string(B))
    elseif typeof(b) <: ChebyshevSpace
        error("Override ConversionOperator from " * string(A) * " to ChebyshevSpace")
    else
        ConversionOperator(a,ChebyshevSpace(AnyDomain()),b)
    end
end

## convert TO Chebyshev
ConversionOperator(A::FunctionSpace)=ConversionOperator(A,ChebyshevSpace(domain(A)))