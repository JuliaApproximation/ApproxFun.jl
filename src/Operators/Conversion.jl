export Conversion

abstract AbstractConversion{T}<:BandedOperator{T}

immutable Conversion{S<:FunctionSpace,V<:FunctionSpace,T} <: AbstractConversion{T}
    domainspace::S
    rangespace::V
end
# Conversion{S<:PeriodicFunctionSpace,V<:PeriodicFunctionSpace}(A::S,B::V)=Conversion{S,V,Complex{Float64}}(A,B)
# Conversion{S<:IntervalFunctionSpace,V<:IntervalFunctionSpace}(A::S,B::V)=Conversion{S,V,Float64}(A,B)

Base.convert{OT<:Operator,S,V}(::Type{OT},C::Conversion{S,V})=Conversion{S,V,eltype(OT)}(C.domainspace,C.rangespace)

domainspace(C::Conversion)=C.domainspace
rangespace(C::Conversion)=C.rangespace




function Conversion(a::FunctionSpace,b::FunctionSpace)
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



## convert TO canonical
Conversion(A::FunctionSpace)=Conversion(A,canonicalspace(A))



## Wrapper
# this allows for a Derivative implementation to return another operator, use a SpaceOperator containing
# the domain and range space
# but continue to know its a derivative

immutable ConversionWrapper{S<:BandedOperator,T} <: AbstractConversion{T}
    op::S
end


ConversionWrapper(B::BandedOperator)=ConversionWrapper{typeof(B),eltype(B)}(B)
Conversion(A::FunctionSpace,B::FunctionSpace,C::FunctionSpace)=ConversionWrapper(TimesOperator(Conversion(B,C),Conversion(A,B)))

# Base.convert{S,T}(::Type{ConversionWrapper{S,T}},D::ConversionWrapper)=ConversionWrapper{S,T}(convert(S,D.op))
Base.convert{CW<:ConversionWrapper}(::Type{CW},D::CW)=D
function Base.convert{BT<:Operator}(::Type{BT},D::ConversionWrapper)
    BO=convert(BandedOperator{eltype(BT)},D.op)
    ConversionWrapper{typeof(BO),eltype(BT)}(BO)::BT
end

addentries!(D::ConversionWrapper,A,k::Range)=addentries!(D.op,A,k)
for func in (:rangespace,:domainspace,:bandinds,:(Base.stride))
    @eval $func(D::ConversionWrapper)=$func(D.op)
end


#TODO: decide
#promotedomainspace(P::Conversion,sp::FunctionSpace)=SpaceOperator(ConstantOperator(one(eltype(P))),sp,sp)
