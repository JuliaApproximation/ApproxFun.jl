export MappedSpace

##Mapped spaces

#TODO: Mapped Fourier
type MappedSpace{S<:FunctionSpace,D<:Domain,T<:Number,DS<:Domain} <: FunctionSpace{T,DS}
    domain::D
    space::S
    MappedSpace(d::D,sp::S)=new(d,sp)
    MappedSpace(d::D)=new(d,S())
    MappedSpace()=new(D(),S())
end

MappedSpace{D<:Domain,T<:Number,DS<:Domain}(d::D,s::FunctionSpace{T,DS})=MappedSpace{typeof(s),D,T,DS}(d,s)

typealias IntervalMappedSpace{S,D} MappedSpace{S,D,Float64,Interval}
typealias PeriodicMappedSpace{S,D,T} MappedSpace{S,D,T,PeriodicInterval}

typealias LineSpace IntervalMappedSpace{Chebyshev,Line}
typealias RaySpace IntervalMappedSpace{Chebyshev,Ray}
typealias CurveSpace{S,T,DS} MappedSpace{S,Curve{S},T,DS}
typealias OpenCurveSpace{S} CurveSpace{S,Float64,Interval}
typealias ClosedCurveSpace{S,T} CurveSpace{S,T,PeriodicInterval}

Space(d::Line)=LineSpace(d)
Space(d::Ray)=RaySpace(d)
Space{S<:FunctionSpace}(d::Curve{S})=CurveSpace{S}(d)


domain(S::MappedSpace)=S.domain
canonicaldomain{D,S}(::Type{IntervalMappedSpace{S,D}})=Interval()
canonicaldomain{D,S,T,DS}(::Type{MappedSpace{S,D,T,DS}})=D()
canonicalspace(S::MappedSpace)=MappedSpace(S.domain,canonicalspace(S.space))

## Construction

Base.ones{T<:Number}(::Type{T},S::MappedSpace)=Fun(ones(T,S.space).coefficients,S)
transform(S::MappedSpace,vals::Vector)=transform(S.space,vals)
itransform(S::MappedSpace,cfs::Vector)=itransform(S.space,cfs)
evaluate{SS,DD,T,TT,DDS}(f::Fun{MappedSpace{SS,DD,TT,DDS},T},x)=evaluate(Fun(coefficients(f),space(f).space),tocanonical(f,x))


for op in (:(Base.first),:(Base.last))
    @eval $op{S<:MappedSpace}(f::Fun{S})=$op(Fun(coefficients(f),space(f).space))
end    



# Transform form chebyshev U series to dirichlet-neumann U series
function uneumann_dirichlet_transform{T<:Number}(v::Vector{T})
    n=length(v)
    w=Array(T,n-4)

    for k = n-4:-1:1
        sc=(3+k)*(4+k)/((-2-k)*(1+k))
        w[k]=sc*v[k+4] 
        
        if k <= n-6
            w[k]-=sc*2*(4+k)/(5+k)*w[k+2]
        end
        if k <= n-8
            w[k]+=sc*((6+k)/(4+k))*w[k+4]
        end
    end
    
    w
end


# This takes a vector in dirichlet-neumann series on [-1,1]
# and return coefficients in T series that satisfy
# (1-x^2)^2 u' = f
function uneumannrange_xsqd{T<:Number}(v::Vector{T})
    n = length(v)
    w=Array(T,n+1)
    
    for k=n:-1:1
        sc=-((16*(1+k)*(2+k))/(k*(3+k)*(4+k)))
        w[k+1]=sc*v[k]
        
        if k <= n-2
            w[k+1]-=sc*(k*(4+k))/(8(k+1))*w[k+3]
        end
        
        if k <= n-4
            w[k+1]+=sc*((k*(k+4))/(16(k+2)))*w[k+5]
        end
    end
    w[1]=zero(T)
    
    w
end




#integration functions
#integration is done by solving (1-x^2)^2 u' = (1-x^2)^2 M' f, u[-1] == 0



function integrate(f::Fun{LineSpace})
    d=domain(f)
    @assert d.α==d.β==-1.
    # || d.α==d.β==-.5
    
#    if domain(f).α==domain(f).β==-1.
        Fun(uneumannrange_xsqd(uneumann_dirichlet_transform(coefficients(Fun([1.5,0.,.5]).*Fun(f.coefficients),UltrasphericalSpace{1}))),f.space)
#    end
#     elseif d.α==d.β==-.5
#         u=divide_singularity(f)
#             integrate(SingFun(Fun(u),-.5,-.5))
#     end  

end

function integrate(f::Fun{RaySpace})
    x=Fun(identity)
    g=fromcanonicalD(f,x)*Fun(f.coefficients)
    Fun(integrate(Fun(g,Chebyshev)).coefficients,space(f))
end

for T in (Float64,Complex{Float64})
    function Base.sum(f::Fun{LineSpace})
        d=domain(f)
        if d.α==d.β==-.5
            sum(Fun(divide_singularity(f.coefficients),JacobiWeightSpace(-.5,-.5,Interval())))
        else
            cf = integrate(f)
            last(cf) - first(cf)
        end
    end
end




## identity

function identity_fun{SS,DD,DDS,DDT}(S::MappedSpace{SS,DD,DDT,DDS})
    sf=fromcanonical(S,Fun(identity,S.space))
    if isa(space(sf),JacobiWeightSpace)
        Fun(coefficients(sf),JacobiWeightSpace(sf.space.α,sf.space.β,S))
    else
         @assert isa(space(sf),S.space)
         Fun(coefficients(sf),S)
    end
end



## Operators

function Evaluation(S1::MappedSpace,x::Bool,order::Integer)
    @assert order==0
    EvaluationWrapper(S1,x,order,Evaluation(S1.space,x,order))
end

Conversion(S1::MappedSpace,S2::MappedSpace)=ConversionWrapper(
    SpaceOperator(Conversion(S1.space,S2.space),
        S1,S2))
        

        
function conversion_rule(S1::MappedSpace,S2::MappedSpace)
    cr=conversion_rule(S1.space,S2.space)
    if cr==S1.space
        S1
    elseif cr==S2.space
        S2
    else
        NoSpace()
    end
end

function addentries!{S1<:MappedSpace,S2<:MappedSpace}(M::Multiplication{S1,S2},A::ShiftArray,kr::Range)
    @assert domain(M.f)==domain(M.space)
    mf=Fun(coefficients(M.f),space(M.f).space)
    addentries!(Multiplication(mf,M.space.space),A,kr)
end



function Derivative(S::MappedSpace,order::Int)
    x=Fun(identity,S)
    D1=Derivative(S.space)
    DS=SpaceOperator(D1,S,MappedSpace(domain(S),rangespace(D1)))
    M=Multiplication(Fun(tocanonicalD(S,x),S),DS|>rangespace)
    D=DerivativeWrapper(M*DS,1)
    if order==1
        D
    else
        Derivative(rangespace(D),order-1)*D
    end
end



## CurveSpace

function evaluate{C<:CurveSpace,T}(f::Fun{C,T},x::Number)
    c=f.space
    rts=roots(domain(f).curve-x)
    @assert length(rts)==1
    evaluate(Fun(f.coefficients,c.space),first(rts))
end


identity_fun{S}(d::CurveSpace{S})=Fun(d.domain.curve.coefficients,d)

