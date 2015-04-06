export Evaluation,ivp

## Evaluation constructors

abstract AbstractEvaluation{T}<:Functional{T}

# M = Bool if endpoint
immutable Evaluation{S<:FunctionSpace,M<:Union(Number,Bool),T<:Number} <: AbstractEvaluation{T}
    space::S
    x::M
    order::Int
end
Evaluation{T}(::Type{T},sp::FunctionSpace,x,order::Integer)=Evaluation{typeof(sp),typeof(x),T}(sp,x,order)
Evaluation(sp::AnySpace,x::Bool,k::Integer)=Evaluation{AnySpace,Bool,UnsetNumber}(sp,x,k)
Evaluation(sp::FunctionSpace{ComplexBasis},x,order::Integer)=Evaluation{typeof(sp),typeof(x),Complex{real(eltype(domain(sp)))}}(sp,x,order)
Evaluation(sp::FunctionSpace,x,order::Integer)=Evaluation{typeof(sp),typeof(x),eltype(domain(sp))}(sp,x,order)

#Evaluation(sp::AnySpace,x::Bool)=Evaluation(sp,x,0)
Evaluation(d::FunctionSpace,x::Union(Number,Bool))=Evaluation(d,x,0)

Evaluation(d::Domain,x::Union(Number,Bool),n...)=Evaluation(Space(d),x,n...)
Evaluation(x::Union(Number,Bool))=Evaluation(AnySpace(),x,0)
Evaluation(x::Union(Number,Bool),k::Integer)=Evaluation(AnySpace(),x,k)
Evaluation{T<:Number}(d::Vector{T},x::Union(Number,Bool),o::Integer)=Evaluation(Interval(d),x,o)


Base.convert{BT<:Operator}(::Type{BT},E::Evaluation)=Evaluation(eltype(BT),E.space,E.x,E.order)


## default getindex
getindex{S,M,T}(D::Evaluation{S,M,T},kr::Range)=T[differentiate(Fun([zeros(T,k-1);one(T)],D.space),D.order)[D.x] for k=kr]

function getindex{S,T}(D::Evaluation{S,Bool,T},kr::Range)
    if !D.x
        T[first(differentiate(Fun([zeros(T,k-1),one(T)],D.space),D.order)) for k=kr]
    else
        T[last(differentiate(Fun([zeros(T,k-1),one(T)],D.space),D.order)) for k=kr]
    end
end




## EvaluationWrapper

immutable EvaluationWrapper{S<:FunctionSpace,M<:Union(Number,Bool),FS<:Functional,T<:Number} <: AbstractEvaluation{T}
    space::S
    x::M
    order::Int
    functional::FS
end

EvaluationWrapper(sp::FunctionSpace,x::Union(Number,Bool),order::Integer,func::Functional)=EvaluationWrapper{typeof(sp),typeof(x),typeof(func),eltype(sp)}(sp,x,order,func)
getindex(E::EvaluationWrapper,kr::Range)=getindex(E.functional,kr)

domainspace(E::AbstractEvaluation)=E.space
domain(E::AbstractEvaluation)=domain(E.space)
promotedomainspace{T}(E::AbstractEvaluation{T},sp::FunctionSpace)=Evaluation(promote_type(T,eltype(sp)),sp,E.x,E.order)
Base.stride(E::EvaluationWrapper)=stride(E.functional)

## Convenience routines


evaluate(d::Domain,x)=Evaluation(d,x)
ldirichlet(d)=Evaluation(d,false)
rdirichlet(d)=Evaluation(d,true)
lneumann(d)=Evaluation(d,false,1)
rneumann(d)=Evaluation(d,true,1)


ldiffbc(d,k) = Evaluation(d,false,k)
rdiffbc(d,k) = Evaluation(d,true,k)


ivp(d)=[ldirichlet(d),lneumann(d)]
dirichlet(d)=[ldirichlet(d),rdirichlet(d)]
neumann(d)=[lneumann(d),rneumann(d)]
diffbcs(d,k) = [ldiffbc(d,k),rdiffbc(d,k)]
periodic(d,k) = Functional{eltype(d)}[Evaluation(d,false,i)-Evaluation(d,true,i) for i=0:k]



for op in (:rdirichlet,:ldirichlet,:dirichlet,:lneumann,:rneumann,:neumann,:ivp)
    @eval begin
        $op()=$op(AnySpace())
        $op(::PeriodicDomain)=[]
    end
end

for op in (:ldiffbc,:rdiffbc,:diffbcs,:periodic)
    @eval $op(k::Integer)=$op(AnySpace(),k)
end




immutable Dirichlet{S,T} <: Operator{T}
    space::S
end
Dirichlet(sp::FunctionSpace)=Dirichlet{typeof(sp),BandedMatrix{eltype(sp)}}(sp)
Dirichlet(d::Domain)=Dirichlet(Space(d))

domainspace(S::Dirichlet)=S.space
rangespace(B::Dirichlet)=Space(∂(domain(B)))

