###
# FunctionBandedOperator allows for defining operators via anonyomous functions
###


immutable FunctionBandedOperator{DS,RS,T} <: Operator{T}
    func::Function
    bandinds::Tuple{Int,Int}
    domainspace::DS
    rangespace::RS
end

FunctionBandedOperator(f,bi,ds,rs)=FunctionBandedOperator{typeof(ds),typeof(rs),typeof(f(1,1))}(f,bi,ds,rs)
FunctionBandedOperator(f,bi)=FunctionBandedOperator(f,bi,AnySpace(),AnySpace())


Base.convert{T}(::Type{Operator{T}},F::FunctionBandedOperator)=FunctionBandedOperator{typeof(F.domainspace),
                                                                    typeof(F.rangespace),
                                                                    T}(F.func,F.bandinds,F.domainspace,F.rangespace)

Base.getindex{DS,RS,T}(F::FunctionBandedOperator{DS,RS,T},k::Integer,j::Integer)=convert(T,F.func(k,j))
bandinds(F::FunctionBandedOperator)=F.bandinds
domainspace(F::FunctionBandedOperator)=F.domainspace
rangespace(F::FunctionBandedOperator)=F.rangespace


BandedOperator(f::Function,bi::Tuple{Int,Int},var...)=FunctionBandedOperator(f,bi,var...)
