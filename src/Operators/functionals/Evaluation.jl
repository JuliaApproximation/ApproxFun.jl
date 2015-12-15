export Evaluation,ivp,bvp

## Evaluation constructors

abstract AbstractEvaluation{T}<:Functional{T}

# M = Bool if endpoint
immutable Evaluation{S,M,T} <: AbstractEvaluation{T}
    space::S
    x::M
    order::Int
end
Evaluation{T}(::Type{T},sp::Space,x::Bool,order::Integer)=Evaluation{typeof(sp),typeof(x),T}(sp,x,order)
function Evaluation{T}(::Type{T},sp::Space,x::Number,order::Integer)
    d=domain(sp)
    if isa(d,IntervalDomain) && isapprox(first(d),x)
        Evaluation(T,sp,false,order)
    elseif isa(d,IntervalDomain) && isapprox(last(d),x)
        Evaluation(T,sp,true,order)
    else
        Evaluation{typeof(sp),typeof(x),T}(sp,x,order)
    end
end

Evaluation(sp::AnySpace,x::Bool,k::Integer)=Evaluation{AnySpace,Bool,UnsetNumber}(sp,x,k)
Evaluation(sp::Space{ComplexBasis},x,order::Integer)=Evaluation(Complex{real(eltype(domain(sp)))},sp,x,order)
Evaluation(sp::Space,x,order::Integer)=Evaluation(eltype(domain(sp)),sp,x,order)

#Evaluation(sp::AnySpace,x::Bool)=Evaluation(sp,x,0)
Evaluation(d::Space,x::Union{Number,Bool})=Evaluation(d,x,0)

Evaluation(d::Domain,x::Union{Number,Bool},n...)=Evaluation(Space(d),x,n...)
Evaluation(x::Union{Number,Bool})=Evaluation(AnySpace(),x,0)
Evaluation(x::Union{Number,Bool},k::Integer)=Evaluation(AnySpace(),x,k)
Evaluation{T<:Number}(d::Vector{T},x::Union{Number,Bool},o::Integer)=Evaluation(Interval(d),x,o)

rangespace{S<:AmbiguousSpace}(E::Evaluation{S,Bool})=ConstantSpace()
rangespace{S}(E::Evaluation{S,Bool})=ConstantSpace(Point(E.x?last(domain(E)):first(domain(E))))
rangespace(E::Evaluation)=ConstantSpace(Point(E.x))

for TYP in (:Operator,:Functional)
    @eval Base.convert{T}(::Type{$TYP{T}},E::Evaluation)=Evaluation(T,E.space,E.x,E.order)
end


## default getindex
getindex{S,M,T}(D::Evaluation{S,M,T},kr::Range)=T[differentiate(Fun([zeros(T,k-1);one(T)],D.space),D.order)(D.x) for k=kr]

function getindex{S,T}(D::Evaluation{S,Bool,T},kr::Range)
    if !D.x
        T[first(differentiate(Fun([zeros(T,k-1);one(T)],D.space),D.order)) for k=kr]
    else
        T[last(differentiate(Fun([zeros(T,k-1);one(T)],D.space),D.order)) for k=kr]
    end
end




## EvaluationWrapper

immutable EvaluationWrapper{S<:Space,M<:Union{Number,Bool},FS<:Functional,T<:Number} <: AbstractEvaluation{T}
    space::S
    x::M
    order::Int
    functional::FS
end

EvaluationWrapper(sp::Space,x::Union{Number,Bool},order::Integer,func::Functional)=EvaluationWrapper{typeof(sp),typeof(x),typeof(func),eltype(sp)}(sp,x,order,func)
getindex(E::EvaluationWrapper,kr::Range)=getindex(E.functional,kr)

domainspace(E::AbstractEvaluation)=E.space
domain(E::AbstractEvaluation)=domain(E.space)
promotedomainspace{T}(E::AbstractEvaluation{T},sp::Space)=Evaluation(promote_type(T,eltype(sp)),sp,E.x,E.order)
Base.stride(E::EvaluationWrapper)=stride(E.functional)

## Convenience routines


evaluate(d::Domain,x)=Evaluation(d,x)
ldiffbc(d,k) = Evaluation(d,false,k)
rdiffbc(d,k) = Evaluation(d,true,k)
diffbcs(d,k) = [ldiffbc(d,k),rdiffbc(d,k)]

ldirichlet(d)=ldiffbc(d,0)
rdirichlet(d)=rdiffbc(d,0)
lneumann(d)=ldiffbc(d,1)
rneumann(d)=rdiffbc(d,1)

dirichlet(d)=diffbcs(d,0)
neumann(d)=diffbcs(d,1)

ivp(d,k) = Functional{eltype(d)}[ldiffbc(d,i) for i=0:k-1]
bvp(d,k) = vcat(Functional{eltype(d)}[ldiffbc(d,i) for i=0:div(k,2)-1],Functional{eltype(d)}[rdiffbc(d,i) for i=0:div(k,2)-1])

periodic(d,k) = Functional{eltype(d)}[Evaluation(d,false,i)-Evaluation(d,true,i) for i=0:k]



for op in (:rdirichlet,:ldirichlet,:dirichlet,:lneumann,:rneumann,:neumann,:ivp,:bvp)
    @eval begin
        $op()=$op(AnySpace())
        $op(::PeriodicDomain)=[]
    end
end

for op in (:ldiffbc,:rdiffbc,:diffbcs,:ivp,:bvp,:periodic)
    @eval $op(k::Integer)=$op(AnySpace(),k)
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
