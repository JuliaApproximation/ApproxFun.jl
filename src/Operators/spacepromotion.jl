

immutable SpaceFunctional{T<:Number,O<:Functional,S<:FunctionSpace} <: Functional{T}
    op::O
    space::S
end

SpaceFunctional{T<:Number,S<:FunctionSpace}(o::Functional{T},s::S)=SpaceFunctional{T,typeof(o),S}(o,s)

getindex(S::SpaceFunctional,k::Range)=getindex(S.op,k)

domainspace(S::SpaceFunctional)=S.space
domain(S::SpaceFunctional)=domain(S.space)

## Space Operator is used to wrap an AnySpace() operator 
immutable SpaceOperator{T<:Number,O<:Operator,S<:FunctionSpace,V<:FunctionSpace} <: BandedOperator{T}
    op::O
    domainspace::S
    rangespace::V
#     
#     function SpaceOperator{T,O,S}(o::O,s::S)
#         @assert domainspace(o)==rangespace(o)==AnySpace()
#         new(o,s)
#     end
end


SpaceOperator{T<:Number}(o::Operator{T},s::FunctionSpace,rs::FunctionSpace)=SpaceOperator{promote_type(T,eltype(s),eltype(rs)),typeof(o),typeof(s),typeof(rs)}(o,s,rs)
SpaceOperator(o,s)=SpaceOperator(o,s,s)

Base.convert{T}(::Type{BandedOperator{T}},S::SpaceOperator)=SpaceOperator(convert(BandedOperator{T},S.op),S.domainspace,S.rangespace)

domain(S::SpaceOperator)=domain(domainspace(S))

domainspace(S::SpaceOperator)=S.domainspace
rangespace(S::SpaceOperator)=S.rangespace
addentries!(S::SpaceOperator,A,kr)=addentries!(S.op,A,kr)

for op in (:bandinds,:(Base.stride))
    @eval $op(S::SpaceOperator)=$op(S.op)
end



##TODO: Do we need both max and min?
function findmindomainspace(ops::Vector)
    sp = AnySpace()
    
    for op in ops
        sp = conversion_type(sp,domainspace(op))
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
        ($op)(P::BandedOperator,::UnsetSpace)=P        
        ($op)(P::BandedOperator,sp::FunctionSpace,::AnySpace)=SpaceOperator(P,sp)
    end
end

promoterangespace(P::Operator,sp::FunctionSpace)=promoterangespace(P,sp,rangespace(P))
promotedomainspace(P::Operator,sp::FunctionSpace)=promotedomainspace(P,sp,domainspace(P))
        
        
promoterangespace(P::Functional,::ConstantSpace,::ConstantSpace)=P # functionals always map to vector space
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

function promotedomainspace{T<:Operator}(ops::Vector{T},S::FunctionSpace)
    k=conversion_type(findmindomainspace(ops),S)
    Operator[promotedomainspace(op,k) for op in ops]
end



####
# choosedomainspace returns a potental domainspace
# where the second argument is a target rangespace
# it defaults to the true domainspace, but if this is ambiguous
# it tries to decide a space.
###

choosedomainspace(A::Operator,::)=domainspace(A)
choosedomainspace(A)=choosedomainspace(A,AnySpace())

function choosedomainspace(ops::Vector,spin)
    sp = AnySpace()
    
    for op in ops
        sp = conversion_type(sp,choosedomainspace(op,spin))
    end
    
    sp
end 


#It's important that domain space is promoted first as it might impact range space
promotespaces(ops::Vector)=promoterangespace(promotedomainspace(ops))
function promotespaces(ops::Vector,b::Fun)
    A=promotespaces(ops)
    if isa(rangespace(A),AmbiguousSpace)
        # try setting the domain space
        A=promoterangespace(promotedomainspace(ops,space(b)))
    end
    A,Fun(b,rangespace(A[end]))
end


