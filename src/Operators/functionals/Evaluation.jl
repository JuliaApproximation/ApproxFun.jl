export Evaluation,ivp,bvp

## Evaluation constructors

abstract Evaluation{T}<:Operator{T}

@functional Evaluation

# M = Bool if endpoint
immutable ConcreteEvaluation{S,M,OT,T} <: Evaluation{T}
    space::S
    x::M
    order::OT
end
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
    differentiate(Fun([zeros(eltype(D),k-1);one(eltype(D))],D.space),D.order)(D.x)


function getindex{S}(D::ConcreteEvaluation{S,Bool},k::Integer)
    T=eltype(D)
    if !D.x
        first(differentiate(Fun([zeros(T,k-1);one(T)],D.space),D.order))
    else
        last(differentiate(Fun([zeros(T,k-1);one(T)],D.space),D.order))
    end
end




## EvaluationWrapper

immutable EvaluationWrapper{S<:Space,M<:Union{Number,Bool},FS<:Operator,OT,T<:Number} <: Evaluation{T}
    space::S
    x::M
    order::OT
    functional::FS
end

EvaluationWrapper(sp::Space,x::Union{Number,Bool},order,func::Operator) =
    EvaluationWrapper{typeof(sp),typeof(x),typeof(func),eltype(sp)}(sp,x,order,func)
getindex(E::EvaluationWrapper,k) = getindex(E.functional,k)

domainspace(E::Evaluation) = E.space
domain(E::Evaluation) = domain(E.space)
promotedomainspace{T}(E::Evaluation{T},sp::Space) =
    Evaluation(promote_type(T,eltype(sp)),sp,E.x,E.order)
Base.stride(E::EvaluationWrapper)=stride(E.functional)

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




immutable Dirichlet{S,T} <: Operator{T}
    space::S
    order::Int
end
Dirichlet(sp::Space)=Dirichlet{typeof(sp),BandedMatrix{eltype(sp)}}(sp,0)
Dirichlet(d::Domain)=Dirichlet(Space(d))
Neumann(sp::Space)=Dirichlet{typeof(sp),BandedMatrix{eltype(sp)}}(sp,1)
Neumann(d::Domain)=Dirichlet(Space(d))


domainspace(S::Dirichlet)=S.space
rangespace(B::Dirichlet)=Space(∂(domain(B)))
