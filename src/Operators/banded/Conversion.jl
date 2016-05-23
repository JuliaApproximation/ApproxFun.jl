export Conversion

abstract Conversion{T}<:BandedOperator{T}

immutable ConcreteConversion{S<:Space,V<:Space,T} <: Conversion{T}
    domainspace::S
    rangespace::V
end

ConcreteConversion(a::Space,b::Space)=ConcreteConversion{typeof(a),typeof(b),
        promote_type(op_eltype_realdomain(a),op_eltype_realdomain(b))}(a,b)


for TYP in (:Operator,:BandedOperator)
    @eval begin
        function Base.convert{T,S,V}(::Type{$TYP{T}},C::ConcreteConversion{S,V})
            if T==eltype(C)
                C
            else
                ConcreteConversion{S,V,T}(C.domainspace,C.rangespace)
            end
        end
    end
end

domainspace(C::ConcreteConversion)=C.domainspace
rangespace(C::ConcreteConversion)=C.rangespace




function defaultConversion(a::Space,b::Space)
    if a==b
        Conversion(a)
    elseif conversion_type(a,b)==NoSpace()
        sp=canonicalspace(a)
        if typeof(sp) == typeof(a)
            error("Implement Conversion from " * string(typeof(sp)) * " to " * string(typeof(b)))
        elseif typeof(sp) == typeof(b)
            error("Implement Conversion from " * string(typeof(a)) * " to " * string(typeof(sp)))
        else
            Conversion(a,sp,b)
        end
    else
        error("Implement Conversion from " * string(typeof(a)) * " to " * string(typeof(b)))
    end
end

Conversion(a::Space,b::Space)=defaultConversion(a,b)
Conversion(a::Space)=ConversionWrapper(eye(a))
Conversion()=ConversionWrapper(eye(UnsetSpace()))


## Wrapper
# this allows for a Derivative implementation to return another operator, use a SpaceOperator containing
# the domain and range space
# but continue to know its a derivative

immutable ConversionWrapper{S<:BandedOperator,T} <: Conversion{T}
    op::S
end

@wrapper ConversionWrapper


ConversionWrapper{T}(::Type{T},op)=ConversionWrapper{typeof(op),T}(op)
ConversionWrapper(B::BandedOperator)=ConversionWrapper{typeof(B),eltype(B)}(B)
Conversion(A::Space,B::Space,C::Space)=Conversion(B,C)*Conversion(A,B)

==(A::ConversionWrapper,B::ConversionWrapper)=A.op==B.op

# Base.convert{S,T}(::Type{ConversionWrapper{S,T}},D::ConversionWrapper)=ConversionWrapper{S,T}(convert(S,D.op))
# Base.convert{CW<:ConversionWrapper}(::Type{CW},D::CW)=D
for TYP in (:Operator,:BandedOperator)
    @eval begin
        function Base.convert{T}(::Type{$TYP{T}},D::ConversionWrapper)
            if T==eltype(D)
                D
            else
                BO=convert(BandedOperator{T},D.op)
                ConversionWrapper{typeof(BO),T}(BO)
            end
        end
    end
end


#promotedomainspace(P::Conversion,sp::Space)=ConversionWrapper(eye(sp))

