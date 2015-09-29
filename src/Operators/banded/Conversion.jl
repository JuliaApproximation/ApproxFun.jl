export Conversion

abstract AbstractConversion{T}<:BandedOperator{T}

immutable Conversion{S<:Space,V<:Space,T} <: AbstractConversion{T}
    domainspace::S
    rangespace::V
end

for TYP in (:Operator,:BandedOperator)
    @eval begin
        function Base.convert{T,S,V}(::Type{$TYP{T}},C::Conversion{S,V})
            if T==eltype(C)
                C
            else
                Conversion{S,V,T}(C.domainspace,C.rangespace)
            end
        end
    end
end

domainspace(C::Conversion)=C.domainspace
rangespace(C::Conversion)=C.rangespace




function defaultconversion(a::Space,b::Space)
    if a==b
        eye(a)
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
        Conversion{typeof(a),typeof(b),promote_type(eltype(a),eltype(b),real(eltype(domain(a))),real(eltype(domain(b))))}(a,b)
    end
end

Conversion(a::Space,b::Space)=defaultconversion(a,b)



## convert TO canonical
Conversion(A::Space)=Conversion(A,canonicalspace(A))



## Wrapper
# this allows for a Derivative implementation to return another operator, use a SpaceOperator containing
# the domain and range space
# but continue to know its a derivative

immutable ConversionWrapper{S<:BandedOperator,T} <: AbstractConversion{T}
    op::S
end


ConversionWrapper(B::BandedOperator)=ConversionWrapper{typeof(B),eltype(B)}(B)
Conversion(A::Space,B::Space,C::Space)=ConversionWrapper(TimesOperator(Conversion(B,C),Conversion(A,B)))

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

addentries!(D::ConversionWrapper,A,k::Range,::Colon)=addentries!(D.op,A,k,:)
for func in (:rangespace,:domainspace,:bandinds,:(Base.stride))
    @eval $func(D::ConversionWrapper)=$func(D.op)
end


#TODO: decide
#promotedomainspace(P::Conversion,sp::Space)=SpaceOperator(ConstantOperator(one(eltype(P))),sp,sp)
