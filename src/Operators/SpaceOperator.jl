




## promotion


## Space Operator is used to wrap an AnySpace() operator 
type SpaceOperator{T<:Number,O<:Operator{T},S<:FunctionSpace} <: BandedOperator{T}
    op::O
    space::S
#     
#     function SpaceOperator{T,O,S}(o::O,s::S)
#         @assert domainspace(o)==rangespace(o)==AnySpace()
#         new(o,s)
#     end
end

SpaceOperator{T<:Number,S<:FunctionSpace}(o::Operator{T},s::S)=SpaceOperator{T,typeof(o),S}(o,s)

domainspace(S::SpaceOperator)=S.space
rangespace(S::SpaceOperator)=S.space
addentries!(S::SpaceOperator,A,kr)=addentries!(S.op,A,kr)
bandinds(S::SpaceOperator)=bandinds(S.op)
domain(S::SpaceOperator)=domain(S.space)


for op in (:promoterangespace,:promotedomainspace)
    @eval begin
        ($op)(P::Operator,::AnySpace)=P
        ($op)(P::Operator,sp::FunctionSpace,::AnySpace)=SpaceOperator(P,sp)
    end
end

promoterangespace(P::Operator,sp::FunctionSpace)=promoterangespace(P,sp,rangespace(P))
promotedomainspace(P::Operator,sp::FunctionSpace)=promotedomainspace(P,sp,domainspace(P))
        
        
promoterangespace(P::Operator,sp::FunctionSpace,cursp::FunctionSpace)=(sp==cursp)?P:TimesOperator(ConversionOperator(cursp,sp),P)
promotedomainspace(P::Operator,sp::FunctionSpace,cursp::FunctionSpace)=(sp==cursp)?P:TimesOperator(P,ConversionOperator(sp,cursp))




function promoterangespace{T<:Operator}(ops::Vector{T})
    k=findmaxrangespace(ops)
    Operator[promoterangespace(op,k) for op in ops]
end




function promotedomainspace{T<:Operator}(ops::Vector{T})
    k=findmindomainspace(ops)
    Operator[promotedomainspace(op,k) for op in ops]
end

#It's important that domain space is promoted first as it might impact range space
promotespaces(ops::Vector)=promoterangespace(promotedomainspace(ops))

