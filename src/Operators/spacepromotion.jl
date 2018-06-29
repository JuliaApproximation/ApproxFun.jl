export →


## Space Operator is used to wrap other operators
# and change the domain/range space
struct SpaceOperator{O<:Operator,S<:Space,V<:Space,T} <: Operator{T}
    op::O
    domainspace::S
    rangespace::V
end

# The promote_type is needed to fix a bug in promotetimes
# not sure if its the right long term solution
SpaceOperator(o::Operator,s::Space,rs::Space) =
    SpaceOperator{typeof(o),typeof(s),typeof(rs),eltype(o)}(o,s,rs)
SpaceOperator(o,s) = SpaceOperator(o,s,s)

function convert(::Type{Operator{T}},S::SpaceOperator) where T
    if T==eltype(S)
        S
    else
        op=convert(Operator{T},S.op)
        SpaceOperator{typeof(op),typeof(S.domainspace),typeof(S.rangespace),T}(op,S.domainspace,S.rangespace)
    end
end



# Similar to wrapper, but different domain/domainspace/rangespace

@wrappergetindex SpaceOperator

# SpaceOperator can change blocks, so we need to override this
getindex(A::SpaceOperator,KR::BlockRange, JR::BlockRange) = defaultgetindex(A,KR,JR)


getindex(A::SpaceOperator,K::Block,J::Block) = A[blockrows(A,K),blockcols(A,J)]
getindex(A::SpaceOperator,K::Block,j) = A[blockrows(A,K),j]
getindex(A::SpaceOperator,k,J::Block) = A[k,blockcols(A,J)]



domain(S::SpaceOperator) = domain(domainspace(S))
domainspace(S::SpaceOperator) = S.domainspace
rangespace(S::SpaceOperator) = S.rangespace



##TODO: Do we need both max and min?
function findmindomainspace(ops::AbstractVector)
    sp = UnsetSpace()

    for op in ops
        sp = union(sp,domainspace(op))
    end

    sp
end

function findmaxrangespace(ops::AbstractVector)
    sp = UnsetSpace()

    for op in ops
        sp = maxspace(sp,rangespace(op))
    end

    sp
end


# The coolest definitions ever!!
# supports Derivative():Chebyshev()→Ultraspherical(1)
Base.colon(A::Operator,b::Space) = promotedomainspace(A,b)
→(A::Operator,b::Space) = promoterangespace(A,b)
Base.colon(A::UniformScaling,b::Space) = Operator(A) : b
→(A::UniformScaling,b::Space) = Operator(A) → b


promoterangespace(P::Operator,sp::Space) = promoterangespace(P,sp,rangespace(P))
promotedomainspace(P::Operator,sp::Space) = promotedomainspace(P,sp,domainspace(P))


promoterangespace(P::Operator,sp::Space,cursp::Space) =
    (sp==cursp) ? P : Conversion(cursp,sp)*P
promotedomainspace(P::Operator,sp::Space,cursp::Space) =
    (sp==cursp) ? P : P*Conversion(sp,cursp)





function promoterangespace(ops::AbstractVector{O}) where O<:Operator
    isempty(ops) && return ops
    k=findmaxrangespace(ops)
    #TODO: T might be incorrect
    T=mapreduce(eltype,promote_type,ops)
    Operator{T}[promoterangespace(op,k) for op in ops]
end
function promotedomainspace(ops::AbstractVector{O}) where O<:Operator
    isempty(ops) && return ops
    k=findmindomainspace(ops)
    #TODO: T might be incorrect
    T=mapreduce(eltype,promote_type,ops)
    Operator{T}[promotedomainspace(op,k) for op in ops]
end
function promotedomainspace(ops::AbstractVector{O},S::Space) where O<:Operator
    isempty(ops) && return ops
    k=conversion_type(findmindomainspace(ops),S)
    #TODO: T might be incorrect
    T=promote_type(mapreduce(eltype,promote_type,ops),prectype(S))
    Operator{T}[promotedomainspace(op,k) for op in ops]
end



####
# choosedomainspace returns a potental domainspace
# where the second argument is a target rangespace
# it defaults to the true domainspace, but if this is ambiguous
# it tries to decide a space.
###

function default_choosedomainspace(A::Operator,sp::Space)
    sp2=domainspace(A)
    isambiguous(sp2) ? sp : sp2
end

choosedomainspace(A::Operator,sp::Space) = default_choosedomainspace(A,sp)

choosedomainspace(A::Operator,f::Fun) = choosedomainspace(A,space(f))
choosedomainspace(A::Operator,f::AbstractVector{FF}) where {FF<:Fun} =
    choosedomainspace(A,Fun(f))
choosedomainspace(A::Operator,_) = choosedomainspace(A)

choosedomainspace(A) = choosedomainspace(A,UnsetSpace())

function choosedomainspace(ops::AbstractVector,spin)
    sp = UnsetSpace()

    for op in ops
        sp = conversion_type(sp,choosedomainspace(op,spin))
    end

    sp
end

choosespaces(A::Operator,b) = promotedomainspace(A,choosedomainspace(A,b))


spacescompatible(A::Operator,B::Operator) =
    spacescompatible(domainspace(A),domainspace(B)) &&
    spacescompatible(rangespace(A),rangespace(B))


#It's important that domain space is promoted first as it might impact range space
promotespaces(ops::AbstractVector) = promoterangespace(promotedomainspace(ops))
function promotespaces(ops::AbstractVector,b::Fun)
    A=promotespaces(ops)
    if isa(rangespace(A),AmbiguousSpace)
        # try setting the domain space
        A=promoterangespace(promotedomainspace(ops,space(b)))
    end
    A,Fun(b,rangespace(A[end]))
end


function promotespaces(A::Operator,B::Operator)
    if spacescompatible(A,B)
        A,B
    else
        tuple(promotespaces([A,B])...)
    end
end




## algebra


A_ldiv_B_coefficients(A::SpaceOperator,b;kwds...) =
    A_ldiv_B_coefficients(A.op,b;kwds...)
