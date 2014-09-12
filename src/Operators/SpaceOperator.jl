

immutable SpaceFunctional{T<:Number,O<:Functional{T},S<:FunctionSpace} <: Functional{T}
    op::O
    space::S
end

SpaceFunctional{T<:Number,S<:FunctionSpace}(o::Functional{T},s::S)=SpaceFunctional{T,typeof(o),S}(o,s)

getindex(S::SpaceFunctional,k::Range)=getindex(S.op,k)


## Space Operator is used to wrap an AnySpace() operator 
immutable SpaceOperator{T<:Number,O<:Operator{T},S<:FunctionSpace} <: BandedOperator{T}
    op::O
    space::S
#     
#     function SpaceOperator{T,O,S}(o::O,s::S)
#         @assert domainspace(o)==rangespace(o)==AnySpace()
#         new(o,s)
#     end
end

SpaceOperator{T<:Number,S<:FunctionSpace}(o::Operator{T},s::S)=SpaceOperator{T,typeof(o),S}(o,s)

for TT in (:SpaceFunctional,:SpaceOperator)
    @eval begin
        domainspace(S::($TT))=S.space
        domain(S::($TT))=domain(S.space)
    end
end


rangespace(S::SpaceOperator)=S.space
addentries!(S::SpaceOperator,A,kr)=addentries!(S.op,A,kr)

bandinds(S::SpaceOperator)=bandinds(S.op)


promotedomainspace(P::Functional,sp::FunctionSpace,::AnySpace)=SpaceFunctional(P,sp)

for op in (:promoterangespace,:promotedomainspace)
    @eval begin
        ($op)(P::Operator,::AnySpace)=P
        ($op)(P::Operator,sp::FunctionSpace,::AnySpace)=SpaceOperator(P,sp)
    end
end

promoterangespace(P::Operator,sp::FunctionSpace)=promoterangespace(P,sp,rangespace(P))
promotedomainspace(P::Operator,sp::FunctionSpace)=promotedomainspace(P,sp,domainspace(P))
        
        

promoterangespace(P::Operator,sp::FunctionSpace,cursp::FunctionSpace)=(sp==cursp)?P:TimesOperator(ConversionOperator(cursp,sp),P)

promotedomainspace(P::Functional,sp::FunctionSpace,cursp::FunctionSpace)=(sp==cursp)?P:TimesFunctional(P,ConversionOperator(sp,cursp))
promotedomainspace(P::Operator,sp::FunctionSpace,cursp::FunctionSpace)=(sp==cursp)?P:TimesOperator(P,ConversionOperator(sp,cursp))





function promoterangespace{T<:Operator}(ops::Vector{T})
    k=findmaxrangespace(ops)
    Operator[promoterangespace(op,k) for op in ops]
end


function promotedomainspace{T<:Functional}(ops::Vector{T})
    k=findmindomainspace(ops)
    Functional[promotedomainspace(op,k) for op in ops]
end


function promotedomainspace{T<:Operator}(ops::Vector{T})
    k=findmindomainspace(ops)
    Operator[promotedomainspace(op,k) for op in ops]
end

#It's important that domain space is promoted first as it might impact range space
promotespaces(ops::Vector)=promoterangespace(promotedomainspace(ops))

