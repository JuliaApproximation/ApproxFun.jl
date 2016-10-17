export Evaluation,ivp,bvp,Dirichlet

## Evaluation constructors

abstract Evaluation{T}<:Operator{T}

@functional Evaluation

# M = Bool if endpoint
immutable ConcreteEvaluation{S,M,OT,T} <: Evaluation{T}
    space::S
    x::M
    order::OT
end

ConcreteEvaluation(sp::Space{RealBasis},x::Number,o::Number) =
    ConcreteEvaluation{typeof(sp),typeof(x),typeof(o),eltype(domain(sp))}(sp,x,o)

Evaluation{T}(::Type{T},sp::UnivariateSpace,x::Bool,order::Integer) =
    ConcreteEvaluation{typeof(sp),typeof(x),typeof(order),T}(sp,x,order)
function Evaluation{T}(::Type{T},sp::UnivariateSpace,x::Number,order::Integer)
    d=domain(sp)
    if isa(d,IntervalDomain) && isapprox(first(d),x)
        Evaluation(T,sp,false,order)
    elseif isa(d,IntervalDomain) && isapprox(last(d),x)
        Evaluation(T,sp,true,order)
    else
        ConcreteEvaluation{typeof(sp),typeof(x),typeof(order),T}(sp,x,order)
    end
end

Evaluation(sp::UnsetSpace,x::Bool,k::Integer) =
    ConcreteEvaluation{UnsetSpace,Bool,typeof(k),UnsetNumber}(sp,x,k)
Evaluation(sp::Space{ComplexBasis},x,order::Integer) =
    Evaluation(Complex{real(eltype(domain(sp)))},sp,x,order)
Evaluation(sp::Space,x,order::Integer) = Evaluation(eltype(domain(sp)),sp,x,order)

#Evaluation(sp::UnsetSpace,x::Bool)=Evaluation(sp,x,0)
Evaluation(d::Space,x::Union{Number,Bool}) = Evaluation(d,x,0)

Evaluation(d::Domain,x::Union{Number,Bool},n...) = Evaluation(Space(d),x,n...)
Evaluation(x::Union{Number,Bool}) = Evaluation(UnsetSpace(),x,0)
Evaluation(x::Union{Number,Bool},k::Integer) = Evaluation(UnsetSpace(),x,k)
Evaluation{T<:Number}(d::Vector{T},x::Union{Number,Bool},o::Integer) = Evaluation(Interval(d),x,o)

rangespace{S<:AmbiguousSpace}(E::ConcreteEvaluation{S,Bool}) = ConstantSpace()
rangespace{S}(E::ConcreteEvaluation{S,Bool}) = ConstantSpace(Point(E.x?last(domain(E)):first(domain(E))))
rangespace(E::ConcreteEvaluation) = ConstantSpace(Point(E.x))


function Base.convert{T}(::Type{Operator{T}},E::ConcreteEvaluation)
    if T == eltype(E)
        E
    else
        ConcreteEvaluation{typeof(E.space),typeof(E.x),typeof(E.order),T}(E.space,E.x,E.order)
    end
end



## default getindex
getindex(D::ConcreteEvaluation,k::Integer) =
    eltype(D)(differentiate(Fun([zeros(eltype(D),k-1);one(eltype(D))],D.space),D.order)(D.x))


function getindex{S}(D::ConcreteEvaluation{S,Bool},k::Integer)
    T=eltype(D)
    if !D.x
        T(first(differentiate(Fun([zeros(T,k-1);one(T)],D.space),D.order)))
    else
        T(last(differentiate(Fun([zeros(T,k-1);one(T)],D.space),D.order)))
    end
end




## EvaluationWrapper

immutable EvaluationWrapper{S<:Space,M,FS<:Operator,OT,T<:Number} <: Evaluation{T}
    space::S
    x::M
    order::OT
    functional::FS
end


#TODO: @wrapper
EvaluationWrapper(sp::Space,x,order,func::Operator) =
    EvaluationWrapper{typeof(sp),typeof(x),typeof(func),typeof(order),eltype(func)}(sp,x,order,func)
getindex(E::EvaluationWrapper,k) = E.functional[k]

domainspace(E::Evaluation) = E.space
domain(E::Evaluation) = domain(E.space)
promotedomainspace(E::Evaluation,sp::Space) = Evaluation(sp,E.x,E.order)
Base.stride(E::EvaluationWrapper)=stride(E.functional)


function Base.convert{T}(::Type{Operator{T}},E::EvaluationWrapper)
    if T == eltype(E)
        E
    else
        EvaluationWrapper(E.space,E.x,E.order,Operator{T}(E.functional))::Operator{T}
    end
end

## Convenience routines


evaluate(d::Domain,x) = Evaluation(d,x)
ldiffbc(d,k) = Evaluation(d,false,k)
rdiffbc(d,k) = Evaluation(d,true,k)
diffbcs(d,k) = [ldiffbc(d,k),rdiffbc(d,k)]

ldirichlet(d) = ldiffbc(d,0)
rdirichlet(d) = rdiffbc(d,0)
lneumann(d) = ldiffbc(d,1)
rneumann(d) = rdiffbc(d,1)

dirichlet(d) = diffbcs(d,0)
neumann(d) = diffbcs(d,1)

ivp(d,k) = Operator{eltype(d)}[ldiffbc(d,i) for i=0:k-1]
bvp(d,k) = vcat(Operator{eltype(d)}[ldiffbc(d,i) for i=0:div(k,2)-1],
                Operator{eltype(d)}[rdiffbc(d,i) for i=0:div(k,2)-1])

periodic(d,k) = Operator{eltype(d)}[Evaluation(d,false,i)-Evaluation(d,true,i) for i=0:k]



for op in (:rdirichlet,:ldirichlet,:dirichlet,:lneumann,:rneumann,:neumann,:ivp,:bvp)
    @eval begin
        $op()=$op(UnsetSpace())
        $op(::PeriodicDomain)=[]
    end
end

for op in (:ldiffbc,:rdiffbc,:diffbcs,:ivp,:bvp,:periodic)
    @eval $op(k::Integer)=$op(UnsetSpace(),k)
end



abstract Dirichlet{S,T} <: Operator{T}


immutable ConcreteDirichlet{S,V,T} <: Dirichlet{S,T}
    domainspace::S
    rangespace::V
    order::Int
end

ConcreteDirichlet(sp::Space,rs::Space,order) =
    ConcreteDirichlet{typeof(sp),typeof(rs),eltype(sp)}(sp,rs,order)
ConcreteDirichlet(sp::Space,order) = ConcreteDirichlet(sp,Space(∂(domain(sp))),order)

immutable DirichletWrapper{S,T} <: Conversion{T}
    op::S
    order::Int
end

@wrapper DirichletWrapper

DirichletWrapper(B::Operator,λ=0) = DirichletWrapper{typeof(B),eltype(B)}(B,λ)


Dirichlet(sp::Space,λ=0) = error("Override getindex for Dirichlet($sp,$λ)")
Dirichlet(d::Domain,λ...) = Dirichlet(Space(d),λ...)
Neumann(sp::Space) = Dirichlet(sp,1)
Neumann(d::Domain) = Neumann(Space(d))


domainspace(B::ConcreteDirichlet) = B.domainspace
rangespace(B::ConcreteDirichlet) = B.rangespace
