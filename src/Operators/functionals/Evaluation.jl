export Evaluation,ivp,bvp,Dirichlet,Neumann

## Evaluation constructors

abstract type Evaluation{T}<:Operator{T} end

@functional Evaluation

# M = first/last if endpoint
struct ConcreteEvaluation{S,M,OT,T} <: Evaluation{T}
    space::S
    x::M
    order::OT
end

ConcreteEvaluation(sp::Space,x,o::Number) =
    ConcreteEvaluation{typeof(sp),typeof(x),typeof(o),rangetype(sp)}(sp,x,o)


Evaluation(::Type{T},sp::Space,x,order) where {T} =
    ConcreteEvaluation{typeof(sp),typeof(x),typeof(order),T}(sp,x,order)
# TODO: This seems like a bad idea: if you are specifying x, just go with the generic version
function Evaluation(::Type{T},sp::UnivariateSpace,x::Number,order) where {T}
    d=domain(sp)
    if isa(d,IntervalDomain) && isapprox(first(d),x)
        Evaluation(T,sp,first,order)
    elseif isa(d,IntervalDomain) && isapprox(last(d),x)
        Evaluation(T,sp,last,order)
    else
        ConcreteEvaluation{typeof(sp),typeof(x),typeof(order),T}(sp,x,order)
    end
end

Evaluation(sp::Space,x,order) = Evaluation(rangetype(sp),sp,x,order)

Evaluation(d::Space,x::Union{Number,typeof(first),typeof(last)}) = Evaluation(d,x,0)
Evaluation(::Type{T},d::Space,n...) where {T} = error("Override Evaluation for $(typeof(d))")
Evaluation(::Type{T},d,n...) where {T} = Evaluation(T,Space(d),n...)
Evaluation(d,n...) = Evaluation(Space(d),n...)
Evaluation(x::Union{Number,typeof(first),typeof(last)}) = Evaluation(UnsetSpace(),x,0)
Evaluation(x::Union{Number,typeof(first),typeof(last)},k::Integer) =
    Evaluation(UnsetSpace(),x,k)

rangespace(E::ConcreteEvaluation{<:AmbiguousSpace}) = ConstantSpace()
rangespace(E::ConcreteEvaluation) = ConstantSpace(Point(E.x))


function convert(::Type{Operator{T}},E::ConcreteEvaluation) where T
    if T == eltype(E)
        E
    else
        ConcreteEvaluation{typeof(E.space),typeof(E.x),typeof(E.order),T}(E.space,E.x,E.order)
    end
end



## default getindex
getindex(D::ConcreteEvaluation,k::Integer) =
    eltype(D)(differentiate(Fun(D.space,[zeros(eltype(D),k-1);one(eltype(D))]),D.order)(D.x))

#special first/last overrides
for OP in (:first,:last)
    @eval begin
        rangespace(E::ConcreteEvaluation{<:AmbiguousSpace,typeof($OP)}) = UnsetSpace()
        function rangespace(E::ConcreteEvaluation{S,typeof($OP)}) where {S}
            d = domain(domainspace(E))
            isambiguous(d) && return ConstantSpace()
            return ConstantSpace(Point($OP(d)))
        end
        function getindex(D::ConcreteEvaluation{S,typeof($OP)},k::Integer) where {S}
            T=eltype(D)

            T($OP(differentiate(Fun(D.space,[zeros(T,k-1);one(T)]),D.order)))
        end
    end
end






## EvaluationWrapper

struct EvaluationWrapper{S<:Space,M,FS<:Operator,OT,T<:Number} <: Evaluation{T}
    space::S
    x::M
    order::OT
    op::FS
end


@wrapper EvaluationWrapper
EvaluationWrapper(sp::Space,x,order,func::Operator) =
    EvaluationWrapper{typeof(sp),typeof(x),typeof(func),typeof(order),eltype(func)}(sp,x,order,func)


domainspace(E::Evaluation) = E.space
domain(E::Evaluation) = domain(E.space)
promotedomainspace(E::Evaluation,sp::Space) = Evaluation(sp,E.x,E.order)



function convert(::Type{Operator{T}},E::EvaluationWrapper) where T
    if T == eltype(E)
        E
    else
        EvaluationWrapper(E.space,E.x,E.order,convert(Operator{T},E.op))::Operator{T}
    end
end

## Convenience routines


evaluate(d::Domain,x) = Evaluation(d,x)
ldiffbc(d,k) = Evaluation(d,first,k)
rdiffbc(d,k) = Evaluation(d,last,k)

ldirichlet(d) = ldiffbc(d,0)
rdirichlet(d) = rdiffbc(d,0)
lneumann(d) = ldiffbc(d,1)
rneumann(d) = rdiffbc(d,1)


ivp(d,k) = Operator{prectype(d)}[ldiffbc(d,i) for i=0:k-1]
bvp(d,k) = vcat(Operator{prectype(d)}[ldiffbc(d,i) for i=0:div(k,2)-1],
                Operator{prectype(d)}[rdiffbc(d,i) for i=0:div(k,2)-1])

periodic(d,k) = Operator{prectype(d)}[Evaluation(d,first,i)-Evaluation(d,last,i) for i=0:k]



for op in (:rdirichlet,:ldirichlet,:lneumann,:rneumann,:ivp,:bvp)
    @eval begin
        $op() = $op(UnsetSpace())
        $op(::PeriodicDomain) = error("Periodic domains do not have boundaries")
    end
end

for op in (:ldiffbc,:rdiffbc,:ivp,:bvp,:periodic)
    @eval $op(k::Integer) = $op(UnsetSpace(),k)
end



abstract type Dirichlet{S,T} <: Operator{T} end


struct ConcreteDirichlet{S,V,T} <: Dirichlet{S,T}
    domainspace::S
    rangespace::V
    order::Int
end

ConcreteDirichlet(sp::Space,rs::Space,order) =
    ConcreteDirichlet{typeof(sp),typeof(rs),rangetype(sp)}(sp,rs,order)
ConcreteDirichlet(sp::Space,order) = ConcreteDirichlet(sp,Space(∂(domain(sp))),order)
ConcreteDirichlet(sp::Space) = ConcreteDirichlet(sp,0)

convert(::Type{Operator{T}},B::ConcreteDirichlet{S,V}) where {S,V,T} =
    ConcreteDirichlet{S,V,T}(B.domainspace,B.rangespace,B.order)


struct DirichletWrapper{S,T} <: Dirichlet{S,T}
    op::S
    order::Int
end

@wrapper DirichletWrapper

DirichletWrapper(B::Operator,λ=0) = DirichletWrapper{typeof(B),eltype(B)}(B,λ)

convert(::Type{Operator{T}},B::DirichletWrapper) where {T} =
    DirichletWrapper(Operator{T}(B.op),B.order)::Operator{T}

# Default is to use diffbca
default_Dirichlet(sp::Space,λ) = DirichletWrapper([ldiffbc(sp,λ);rdiffbc(sp,λ)],λ)
Dirichlet(sp::Space,λ) = default_Dirichlet(sp,λ)
Dirichlet(sp::Space) = Dirichlet(sp,0)
Dirichlet() = Dirichlet(UnsetSpace())

Dirichlet(d::Domain,λ...) = Dirichlet(Space(d),λ...)
Neumann(sp::Space) = Dirichlet(sp,1)
Neumann(d::Domain) = Neumann(Space(d))
Neumann() = Dirichlet(UnsetSpace(),1)

domainspace(B::ConcreteDirichlet) = B.domainspace
rangespace(B::ConcreteDirichlet) = B.rangespace

promotedomainspace(E::Dirichlet,sp::Space) = Dirichlet(sp,E.order)



"""
`Evaluation(sp,x,k)` is the functional associated with evaluating the
`k`-th derivative at a point `x` for the space `sp`.
"""
Evaluation(sp::Space,x,k)

"""
`Evaluation(sp,x)` is the functional associated with evaluating
at a point `x` for the space `sp`.
"""
Evaluation(sp::Space,x)


"""
`Evaluation(x)` is the functional associated with evaluating
at a point `x`.
"""
Evaluation(x)


"""
`Dirichlet(sp,k)` is the operator associated with restricting the
`k`-th derivative on the boundary for the space `sp`.
"""
Dirichlet(sp::Space,k)

"""
`Dirichlet(sp)` is the operator associated with restricting the
 the boundary for the space `sp`.
"""
Dirichlet(sp::Space)


"""
`Dirichlet()` is the operator associated with restricting on the
 the boundary.
"""
Dirichlet()




"""
`Neumann(sp)` is the operator associated with restricting the
normal derivative on the boundary for the space `sp`.
At the moment it is implemented as `Dirichlet(sp,1)`.
"""
Neumann(sp::Space)

"""
`Neumann( is the operator associated with restricting the
normal derivative on the boundary.
"""
Neumann()
