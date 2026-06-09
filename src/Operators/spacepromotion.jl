export ↦


## Space Operator is used to wrap other operators
# and change the domain/range space
immutable SpaceOperator{O<:Operator,S<:Space,V<:Space,T} <: Operator{T}
    op::O
    domainspace::S
    rangespace::V
end

# The promote_type is needed to fix a bug in promotetimes
# not sure if its the right long term solution
SpaceOperator(o::Operator,s::Space,rs::Space) =
    SpaceOperator{typeof(o),typeof(s),typeof(rs),eltype(o)}(o,s,rs)
SpaceOperator(o,s) = SpaceOperator(o,s,s)

function Base.convert{T}(::Type{Operator{T}},S::SpaceOperator)
    if T==eltype(S)
        S
    else
        op=convert(Operator{T},S.op)
        SpaceOperator{typeof(op),typeof(S.domainspace),typeof(S.rangespace),T}(op,S.domainspace,S.rangespace)
    end
end



# Similar to wrapper, but different domain/domainspace/rangespace

@wrappergetindex SpaceOperator


domain(S::SpaceOperator) = domain(domainspace(S))
domainspace(S::SpaceOperator) = S.domainspace
rangespace(S::SpaceOperator) = S.rangespace



##TODO: Do we need both max and min?
function findmindomainspace(ops::Vector)
    sp = UnsetSpace()

    for op in ops
        sp = conversion_type(sp,domainspace(op))
    end

    sp
end

function findmaxrangespace(ops::Vector)
    sp = UnsetSpace()

    for op in ops
        sp = maxspace(sp,rangespace(op))
    end

    sp
end


# The coolest definitions ever!!
# supports Derivative():Chebyshev()↦Ultraspherical(1)
↦(A::Operator,b::Space) = promoterangespace(A,b)
Base.colon(A::Operator,b::Space) = promotedomainspace(A,b)

promoterangespace(P::Operator,sp::Space) = promoterangespace(P,sp,rangespace(P))
promotedomainspace(P::Operator,sp::Space) = promotedomainspace(P,sp,domainspace(P))


promoterangespace(P::Operator,sp::Space,cursp::Space) =
    (sp==cursp)?P:Conversion(cursp,sp)*P
promotedomainspace(P::Operator,sp::Space,cursp::Space) =
    (sp==cursp)?P:P*Conversion(sp,cursp)





function promoterangespace{O<:Operator}(ops::Vector{O})
    isempty(ops) && return ops
    k=findmaxrangespace(ops)
    #TODO: T might be incorrect
    T=mapreduce(eltype,promote_type,ops)
    Operator{T}[promoterangespace(op,k) for op in ops]
end
function promotedomainspace{O<:Operator}(ops::Vector{O})
    isempty(ops) && return ops
    k=findmindomainspace(ops)
    #TODO: T might be incorrect
    T=mapreduce(eltype,promote_type,ops)
    Operator{T}[promotedomainspace(op,k) for op in ops]
end
function promotedomainspace{O<:Operator}(ops::Vector{O},S::Space)
    isempty(ops) && return ops
    k=conversion_type(findmindomainspace(ops),S)
    #TODO: T might be incorrect
    T=promote_type(mapreduce(eltype,promote_type,ops),eltype(S))
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
    isambiguous(sp2)?sp:sp2
end

choosedomainspace(A::Operator,sp::Space) = default_choosedomainspace(A,sp)

choosedomainspace(A::Operator,f::Fun) = choosedomainspace(A,space(f))
choosedomainspace{FF<:Fun}(A::Operator,f::Vector{FF}) =
    choosedomainspace(A,devec(f))
choosedomainspace(A::Operator,::) = choosedomainspace(A)

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


linsolve(A::SpaceOperator,b::Fun;kwds...) =
    setspace(linsolve(A.op,coefficients(b,rangespace(A));kwds...),domainspace(A))

linsolve{T<:Number}(A::SpaceOperator,b::Array{T};kwds...) =
    setspace(linsolve(A.op,b;kwds...),domainspace(A))
linsolve(A::SpaceOperator,b::Number;kwds...) =
    setspace(linsolve(A.op,b;kwds...),domainspace(A))
