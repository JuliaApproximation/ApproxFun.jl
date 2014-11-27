##Line calculus


type LineSpace <: IntervalDomainSpace
    domain::Line
end

type RaySpace <: IntervalDomainSpace
    domain::Ray
end


Space(d::Line)=LineSpace(d)
Space(d::Ray)=RaySpace(d)


canonicaldomain{T<:LineSpace}(::Type{T})=Line()
canonicaldomain{T<:RaySpace}(::Type{T})=Ray()


## Construction

#domain(S) may be any domain
for TYP in (:RaySpace,:LineSpace)
    @eval begin
        Base.ones{T<:Number}(::Type{T},S::$TYP)=Fun(ones(T,1),S)
        transform(::$TYP,vals::Vector)=chebyshevtransform(vals)
        itransform(::$TYP,cfs::Vector)=ichebyshevtransform(cfs)        
        evaluate(f::Fun{$TYP},x)=clenshaw(f.coefficients,tocanonical(f,x))
    end
    
    for op in (:(Base.first),:(Base.last))
        @eval $op(f::Fun{$TYP})=$op(Fun(f.coefficients))
    end    
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
    f=Fun(x->exp(-x^2),[0.,Inf])
    x=Fun(identity)
    g=fromcanonicalD(f,x)*Fun(f.coefficients)
    Fun(integrate(Fun(g,ChebyshevSpace)).coefficients,space(f))
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

function identity_fun(S::Union(LineSpace,RaySpace))
    sf=fromcanonical(S,Fun(identity))
    Fun(coefficients(sf),JacobiWeightSpace(sf.space.α,sf.space.β,S))
end



## Operators

addentries!{S<:Union(RaySpace,LineSpace)}(M::Multiplication{S,S},A::ShiftArray,kr::Range)=chebmult_addentries!(coefficients(M.f),A,kr)

