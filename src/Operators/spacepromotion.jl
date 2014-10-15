

immutable SpaceFunctional{T<:Number,O<:Functional{T},S<:FunctionSpace} <: Functional{T}
    op::O
    space::S
end

SpaceFunctional{T<:Number,S<:FunctionSpace}(o::Functional{T},s::S)=SpaceFunctional{T,typeof(o),S}(o,s)

getindex(S::SpaceFunctional,k::Range)=getindex(S.op,k)

domainspace(S::SpaceFunctional)=S.space
domain(S::SpaceFunctional)=domain(S.space)

## Space Operator is used to wrap an AnySpace() operator 
immutable SpaceOperator{T<:Number,O<:Operator{T},S<:FunctionSpace,V<:FunctionSpace} <: BandedOperator{T}
    op::O
    domainspace::S
    rangespace::V
#     
#     function SpaceOperator{T,O,S}(o::O,s::S)
#         @assert domainspace(o)==rangespace(o)==AnySpace()
#         new(o,s)
#     end
end


SpaceOperator{T<:Number,S<:FunctionSpace,V<:FunctionSpace}(o::Operator{T},s::S,rs::V)=SpaceOperator{T,typeof(o),S,V}(o,s,rs)
SpaceOperator(o,s)=SpaceOperator(o,s,s)


domain(S::SpaceOperator)=domain(domainspace(S))

domainspace(S::SpaceOperator)=S.domainspace
rangespace(S::SpaceOperator)=S.rangespace
addentries!(S::SpaceOperator,A,kr)=addentries!(S.op,A,kr)

bandinds(S::SpaceOperator)=bandinds(S.op)



##TODO: Do we need both max and min?
function findmindomainspace(ops::Vector)
    sp = AnySpace()
    
    for op in ops
        sp = minspace(sp,domainspace(op))
    end
    
    sp
end

function findmaxrangespace(ops::Vector)
    sp = AnySpace()
    
    for op in ops
        sp = maxspace(sp,rangespace(op))
    end
    
    sp
end




promotedomainspace(P::Functional,sp::FunctionSpace,::AnySpace)=SpaceFunctional(P,sp)

for op in (:promoterangespace,:promotedomainspace)
    @eval begin
        ($op)(P::BandedOperator,::AnySpace)=P
        ($op)(P::BandedOperator,sp::FunctionSpace,::AnySpace)=SpaceOperator(P,sp)
    end
end

promoterangespace(P::Operator,sp::FunctionSpace)=promoterangespace(P,sp,rangespace(P))
promotedomainspace(P::Operator,sp::FunctionSpace)=promotedomainspace(P,sp,domainspace(P))
        
        
promoterangespace(P::Functional,::VectorSpace{1},::VectorSpace{1})=P # functionals always map to vector space
promoterangespace(P::BandedOperator,sp::FunctionSpace,cursp::FunctionSpace)=(sp==cursp)?P:TimesOperator(Conversion(cursp,sp),P)
promotedomainspace(P::Functional,sp::FunctionSpace,cursp::FunctionSpace)=(sp==cursp)?P:TimesFunctional(P,Conversion(sp,cursp))
promotedomainspace(P::BandedOperator,sp::FunctionSpace,cursp::FunctionSpace)=(sp==cursp)?P:TimesOperator(P,Conversion(sp,cursp))





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

